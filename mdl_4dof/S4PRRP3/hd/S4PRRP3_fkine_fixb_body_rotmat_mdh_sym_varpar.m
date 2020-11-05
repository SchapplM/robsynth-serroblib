% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:36
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:36:41
	% EndTime: 2020-11-04 19:36:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:36:41
	% EndTime: 2020-11-04 19:36:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(pkin(6));
	t26 = sin(pkin(6));
	t1 = [t27, -t26, 0, 0; t26, t27, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:36:41
	% EndTime: 2020-11-04 19:36:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t30 = pkin(6) + qJ(2);
	t29 = cos(t30);
	t28 = sin(t30);
	t1 = [t29, -t28, 0, cos(pkin(6)) * pkin(1) + 0; t28, t29, 0, sin(pkin(6)) * pkin(1) + 0; 0, 0, 1, pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:36:41
	% EndTime: 2020-11-04 19:36:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t35 = cos(qJ(3));
	t34 = sin(qJ(3));
	t33 = pkin(6) + qJ(2);
	t32 = cos(t33);
	t31 = sin(t33);
	t1 = [t32 * t35, -t32 * t34, t31, t32 * pkin(2) + t31 * pkin(5) + cos(pkin(6)) * pkin(1) + 0; t31 * t35, -t31 * t34, -t32, t31 * pkin(2) - t32 * pkin(5) + sin(pkin(6)) * pkin(1) + 0; t34, t35, 0, pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:36:41
	% EndTime: 2020-11-04 19:36:41
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (26->15), mult. (13->12), div. (0->0), fcn. (21->6), ass. (0->8)
	t42 = cos(qJ(3));
	t41 = sin(qJ(3));
	t40 = -qJ(4) - pkin(5);
	t39 = pkin(6) + qJ(2);
	t38 = cos(t39);
	t37 = sin(t39);
	t36 = t42 * pkin(3) + pkin(2);
	t1 = [t38 * t42, -t38 * t41, t37, t38 * t36 - t37 * t40 + cos(pkin(6)) * pkin(1) + 0; t37 * t42, -t37 * t41, -t38, t37 * t36 + t38 * t40 + sin(pkin(6)) * pkin(1) + 0; t41, t42, 0, t41 * pkin(3) + pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end