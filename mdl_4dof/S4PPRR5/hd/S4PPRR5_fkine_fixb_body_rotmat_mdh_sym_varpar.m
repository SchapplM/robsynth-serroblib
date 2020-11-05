% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PPRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:33
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:08
	% EndTime: 2020-11-04 19:33:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:08
	% EndTime: 2020-11-04 19:33:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t28 = cos(pkin(6));
	t27 = sin(pkin(6));
	t1 = [t28, -t27, 0, 0; t27, t28, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:08
	% EndTime: 2020-11-04 19:33:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (8->8), mult. (4->4), div. (0->0), fcn. (8->2), ass. (0->3)
	t30 = cos(pkin(6));
	t29 = sin(pkin(6));
	t1 = [0, -t30, t29, t30 * pkin(1) + t29 * qJ(2) + 0; 0, -t29, -t30, t29 * pkin(1) - t30 * qJ(2) + 0; 1, 0, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:08
	% EndTime: 2020-11-04 19:33:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->10), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->6)
	t35 = pkin(1) + pkin(4);
	t34 = cos(qJ(3));
	t33 = sin(qJ(3));
	t32 = cos(pkin(6));
	t31 = sin(pkin(6));
	t1 = [t31 * t33, t31 * t34, t32, t31 * qJ(2) + t35 * t32 + 0; -t32 * t33, -t32 * t34, t31, -t32 * qJ(2) + t35 * t31 + 0; t34, -t33, 0, pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:33:08
	% EndTime: 2020-11-04 19:33:08
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (20->17), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->11)
	t38 = sin(qJ(4));
	t39 = sin(qJ(3));
	t45 = t38 * t39;
	t40 = cos(qJ(4));
	t44 = t39 * t40;
	t41 = cos(qJ(3));
	t43 = pkin(3) * t39 - pkin(5) * t41 + qJ(2);
	t42 = pkin(1) + pkin(4);
	t37 = cos(pkin(6));
	t36 = sin(pkin(6));
	t1 = [t36 * t44 + t37 * t38, -t36 * t45 + t37 * t40, -t36 * t41, t43 * t36 + t42 * t37 + 0; t36 * t38 - t37 * t44, t36 * t40 + t37 * t45, t37 * t41, t42 * t36 - t43 * t37 + 0; t41 * t40, -t41 * t38, t39, t41 * pkin(3) + t39 * pkin(5) + pkin(2) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end