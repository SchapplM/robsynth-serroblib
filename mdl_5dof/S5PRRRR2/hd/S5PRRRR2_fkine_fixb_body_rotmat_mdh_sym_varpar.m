% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:07
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [6x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:50
	% EndTime: 2020-11-04 20:07:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:50
	% EndTime: 2020-11-04 20:07:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:50
	% EndTime: 2020-11-04 20:07:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t32 = cos(qJ(2));
	t31 = sin(qJ(2));
	t1 = [t32, -t31, 0, pkin(1) + 0; t31, t32, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:50
	% EndTime: 2020-11-04 20:07:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (10->7), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t35 = qJ(2) + qJ(3);
	t34 = cos(t35);
	t33 = sin(t35);
	t1 = [t34, -t33, 0, cos(qJ(2)) * pkin(2) + pkin(1) + 0; t33, t34, 0, sin(qJ(2)) * pkin(2) + 0; 0, 0, 1, pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:50
	% EndTime: 2020-11-04 20:07:50
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->11), mult. (4->4), div. (0->0), fcn. (8->6), ass. (0->5)
	t39 = qJ(2) + qJ(3);
	t38 = qJ(4) + t39;
	t37 = cos(t38);
	t36 = sin(t38);
	t1 = [t37, -t36, 0, pkin(3) * cos(t39) + cos(qJ(2)) * pkin(2) + pkin(1) + 0; t36, t37, 0, pkin(3) * sin(t39) + sin(qJ(2)) * pkin(2) + 0; 0, 0, 1, pkin(5) + pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:07:50
	% EndTime: 2020-11-04 20:07:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (31->15), mult. (10->10), div. (0->0), fcn. (18->8), ass. (0->7)
	t43 = qJ(2) + qJ(3);
	t45 = cos(qJ(5));
	t44 = sin(qJ(5));
	t42 = qJ(4) + t43;
	t41 = cos(t42);
	t40 = sin(t42);
	t1 = [t41 * t45, -t41 * t44, t40, t40 * pkin(6) + pkin(3) * cos(t43) + cos(qJ(2)) * pkin(2) + pkin(1) + 0; t40 * t45, -t40 * t44, -t41, -t41 * pkin(6) + pkin(3) * sin(t43) + sin(qJ(2)) * pkin(2) + 0; t44, t45, 0, pkin(5) + pkin(4) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end