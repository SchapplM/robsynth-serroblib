% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPRP1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:31
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [5x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:31:13
	% EndTime: 2020-11-04 19:31:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:31:13
	% EndTime: 2020-11-04 19:31:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(pkin(5));
	t26 = sin(pkin(5));
	t1 = [t27, -t26, 0, 0; t26, t27, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:31:13
	% EndTime: 2020-11-04 19:31:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (6->6), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t29 = cos(pkin(5));
	t28 = sin(pkin(5));
	t1 = [t29, 0, t28, t29 * pkin(1) + t28 * qJ(2) + 0; t28, 0, -t29, t28 * pkin(1) - t29 * qJ(2) + 0; 0, 1, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:31:13
	% EndTime: 2020-11-04 19:31:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (12->10), mult. (12->8), div. (0->0), fcn. (20->4), ass. (0->8)
	t36 = pkin(1) + pkin(2);
	t35 = cos(qJ(3));
	t34 = sin(qJ(3));
	t33 = cos(pkin(5));
	t32 = sin(pkin(5));
	t31 = t32 * t35 - t33 * t34;
	t30 = -t32 * t34 - t33 * t35;
	t1 = [-t30, t31, 0, t32 * qJ(2) + t36 * t33 + 0; t31, t30, 0, -t33 * qJ(2) + t36 * t32 + 0; 0, 0, -1, -pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:31:13
	% EndTime: 2020-11-04 19:31:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (20->16), mult. (24->16), div. (0->0), fcn. (32->4), ass. (0->10)
	t45 = pkin(1) + pkin(2);
	t44 = cos(qJ(3));
	t43 = sin(qJ(3));
	t42 = cos(pkin(5));
	t41 = sin(pkin(5));
	t40 = t42 * pkin(3) - qJ(4) * t41;
	t39 = pkin(3) * t41 + t42 * qJ(4);
	t38 = t41 * t43 + t42 * t44;
	t37 = -t41 * t44 + t42 * t43;
	t1 = [t38, 0, t37, t41 * qJ(2) + t39 * t43 + t40 * t44 + t45 * t42 + 0; -t37, 0, t38, -t42 * qJ(2) + t39 * t44 - t40 * t43 + t45 * t41 + 0; 0, -1, 0, -pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end