% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4RRPR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:46
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4RRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [7x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:46:49
	% EndTime: 2020-11-04 19:46:49
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:46:49
	% EndTime: 2020-11-04 19:46:49
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = cos(qJ(1));
	t30 = sin(qJ(1));
	t1 = [t31, -t30, 0, 0; t30, t31, 0, 0; 0, 0, 1, pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:46:49
	% EndTime: 2020-11-04 19:46:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t34 = qJ(1) + qJ(2);
	t33 = cos(t34);
	t32 = sin(t34);
	t1 = [t33, -t32, 0, cos(qJ(1)) * pkin(1) + 0; t32, t33, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, pkin(5) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:46:49
	% EndTime: 2020-11-04 19:46:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (18->10), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t38 = qJ(1) + qJ(2);
	t37 = pkin(7) + t38;
	t36 = cos(t37);
	t35 = sin(t37);
	t1 = [t36, -t35, 0, pkin(2) * cos(t38) + cos(qJ(1)) * pkin(1) + 0; t35, t36, 0, pkin(2) * sin(t38) + sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(3) + pkin(5) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:46:49
	% EndTime: 2020-11-04 19:46:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (36->16), mult. (12->12), div. (0->0), fcn. (20->8), ass. (0->7)
	t42 = qJ(1) + qJ(2);
	t44 = cos(qJ(4));
	t43 = sin(qJ(4));
	t41 = pkin(7) + t42;
	t40 = cos(t41);
	t39 = sin(t41);
	t1 = [t40 * t44, -t40 * t43, t39, t40 * pkin(3) + t39 * pkin(6) + pkin(2) * cos(t42) + cos(qJ(1)) * pkin(1) + 0; t39 * t44, -t39 * t43, -t40, t39 * pkin(3) - t40 * pkin(6) + pkin(2) * sin(t42) + sin(qJ(1)) * pkin(1) + 0; t43, t44, 0, qJ(3) + pkin(5) + pkin(4) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end