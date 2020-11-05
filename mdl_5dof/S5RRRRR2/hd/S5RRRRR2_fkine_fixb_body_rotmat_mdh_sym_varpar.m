% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:48
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [2x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:31
	% EndTime: 2020-11-04 20:48:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:31
	% EndTime: 2020-11-04 20:48:31
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (1->1), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t39 = cos(qJ(1));
	t38 = sin(qJ(1));
	t1 = [t39, -t38, 0, 0; t38, t39, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:31
	% EndTime: 2020-11-04 20:48:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (7->4), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t42 = qJ(1) + qJ(2);
	t41 = cos(t42);
	t40 = sin(t42);
	t1 = [t41, -t40, 0, cos(qJ(1)) * pkin(1) + 0; t40, t41, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:31
	% EndTime: 2020-11-04 20:48:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->6), mult. (6->6), div. (0->0), fcn. (14->6), ass. (0->6)
	t47 = cos(qJ(3));
	t46 = sin(qJ(3));
	t45 = qJ(1) + qJ(2);
	t44 = cos(t45);
	t43 = sin(t45);
	t1 = [t44 * t47, -t44 * t46, t43, cos(qJ(1)) * pkin(1) + 0; t43 * t47, -t43 * t46, -t44, sin(qJ(1)) * pkin(1) + 0; t46, t47, 0, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:31
	% EndTime: 2020-11-04 20:48:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (22->10), mult. (11->10), div. (0->0), fcn. (19->8), ass. (0->8)
	t55 = pkin(2) * cos(qJ(3));
	t53 = qJ(1) + qJ(2);
	t52 = qJ(3) + qJ(4);
	t51 = cos(t53);
	t50 = cos(t52);
	t49 = sin(t53);
	t48 = sin(t52);
	t1 = [t51 * t50, -t51 * t48, t49, t51 * t55 + cos(qJ(1)) * pkin(1) + 0; t49 * t50, -t49 * t48, -t51, t49 * t55 + sin(qJ(1)) * pkin(1) + 0; t48, t50, 0, sin(qJ(3)) * pkin(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:48:31
	% EndTime: 2020-11-04 20:48:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (32->13), mult. (23->18), div. (0->0), fcn. (36->10), ass. (0->14)
	t69 = pkin(2) * cos(qJ(3));
	t61 = qJ(1) + qJ(2);
	t57 = sin(t61);
	t62 = sin(qJ(5));
	t68 = t57 * t62;
	t63 = cos(qJ(5));
	t67 = t57 * t63;
	t59 = cos(t61);
	t66 = t59 * t62;
	t65 = t59 * t63;
	t60 = qJ(3) + qJ(4);
	t58 = cos(t60);
	t56 = sin(t60);
	t1 = [t58 * t65 + t68, -t58 * t66 + t67, t59 * t56, t59 * t69 + cos(qJ(1)) * pkin(1) + 0; t58 * t67 - t66, -t58 * t68 - t65, t57 * t56, t57 * t69 + sin(qJ(1)) * pkin(1) + 0; t56 * t63, -t56 * t62, -t58, sin(qJ(3)) * pkin(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end