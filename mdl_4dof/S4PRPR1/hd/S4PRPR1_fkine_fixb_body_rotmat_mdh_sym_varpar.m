% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRPR1 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:34
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:34:12
	% EndTime: 2020-11-04 19:34:12
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:34:12
	% EndTime: 2020-11-04 19:34:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(pkin(6));
	t26 = sin(pkin(6));
	t1 = [t27, -t26, 0, 0; t26, t27, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:34:13
	% EndTime: 2020-11-04 19:34:13
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
	% StartTime: 2020-11-04 19:34:13
	% EndTime: 2020-11-04 19:34:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (17->10), mult. (6->6), div. (0->0), fcn. (10->4), ass. (0->4)
	t33 = pkin(6) + qJ(2);
	t32 = cos(t33);
	t31 = sin(t33);
	t1 = [t32, 0, t31, t32 * pkin(2) + t31 * qJ(3) + cos(pkin(6)) * pkin(1) + 0; t31, 0, -t32, t31 * pkin(2) - t32 * qJ(3) + sin(pkin(6)) * pkin(1) + 0; 0, 1, 0, pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:34:13
	% EndTime: 2020-11-04 19:34:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (27->14), mult. (14->10), div. (0->0), fcn. (22->6), ass. (0->9)
	t41 = pkin(2) + pkin(3);
	t40 = cos(qJ(4));
	t39 = sin(qJ(4));
	t38 = pkin(6) + qJ(2);
	t37 = cos(t38);
	t36 = sin(t38);
	t35 = t36 * t40 - t37 * t39;
	t34 = -t36 * t39 - t37 * t40;
	t1 = [-t34, t35, 0, t41 * t37 + t36 * qJ(3) + cos(pkin(6)) * pkin(1) + 0; t35, t34, 0, t41 * t36 - t37 * qJ(3) + sin(pkin(6)) * pkin(1) + 0; 0, 0, -1, -pkin(5) + pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end