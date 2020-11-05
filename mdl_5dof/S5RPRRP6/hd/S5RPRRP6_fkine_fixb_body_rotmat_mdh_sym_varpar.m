% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPRRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:23
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:58
	% EndTime: 2020-11-04 20:23:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:58
	% EndTime: 2020-11-04 20:23:58
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [t46, -t45, 0, 0; t45, t46, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:58
	% EndTime: 2020-11-04 20:23:58
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (9->6), mult. (2->2), div. (0->0), fcn. (6->4), ass. (0->4)
	t49 = qJ(1) + pkin(8);
	t48 = cos(t49);
	t47 = sin(t49);
	t1 = [t48, -t47, 0, cos(qJ(1)) * pkin(1) + 0; t47, t48, 0, sin(qJ(1)) * pkin(1) + 0; 0, 0, 1, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:58
	% EndTime: 2020-11-04 20:23:58
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (21->12), mult. (10->10), div. (0->0), fcn. (18->6), ass. (0->6)
	t54 = cos(qJ(3));
	t53 = sin(qJ(3));
	t52 = qJ(1) + pkin(8);
	t51 = cos(t52);
	t50 = sin(t52);
	t1 = [t51 * t54, -t51 * t53, t50, t51 * pkin(2) + t50 * pkin(6) + cos(qJ(1)) * pkin(1) + 0; t50 * t54, -t50 * t53, -t51, t50 * pkin(2) - t51 * pkin(6) + sin(qJ(1)) * pkin(1) + 0; t53, t54, 0, qJ(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:58
	% EndTime: 2020-11-04 20:23:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (38->19), mult. (32->24), div. (0->0), fcn. (45->8), ass. (0->11)
	t58 = sin(qJ(4));
	t61 = cos(qJ(3));
	t64 = t58 * t61;
	t60 = cos(qJ(4));
	t63 = t60 * t61;
	t59 = sin(qJ(3));
	t62 = pkin(3) * t61 + pkin(7) * t59 + pkin(2);
	t57 = qJ(1) + pkin(8);
	t56 = cos(t57);
	t55 = sin(t57);
	t1 = [t55 * t58 + t56 * t63, t55 * t60 - t56 * t64, t56 * t59, cos(qJ(1)) * pkin(1) + t55 * pkin(6) + 0 + t62 * t56; t55 * t63 - t56 * t58, -t55 * t64 - t56 * t60, t55 * t59, sin(qJ(1)) * pkin(1) - t56 * pkin(6) + 0 + t62 * t55; t59 * t60, -t59 * t58, -t61, t59 * pkin(3) - t61 * pkin(7) + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:23:58
	% EndTime: 2020-11-04 20:23:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (48->22), mult. (39->26), div. (0->0), fcn. (52->8), ass. (0->14)
	t70 = sin(qJ(4));
	t73 = cos(qJ(3));
	t77 = t70 * t73;
	t72 = cos(qJ(4));
	t76 = t72 * t73;
	t75 = pkin(4) * t70 + pkin(6);
	t65 = t72 * pkin(4) + pkin(3);
	t69 = -qJ(5) - pkin(7);
	t71 = sin(qJ(3));
	t74 = t65 * t73 - t69 * t71 + pkin(2);
	t68 = qJ(1) + pkin(8);
	t67 = cos(t68);
	t66 = sin(t68);
	t1 = [t66 * t70 + t67 * t76, t66 * t72 - t67 * t77, t67 * t71, cos(qJ(1)) * pkin(1) + 0 + t75 * t66 + t74 * t67; t66 * t76 - t67 * t70, -t66 * t77 - t67 * t72, t66 * t71, sin(qJ(1)) * pkin(1) + 0 - t75 * t67 + t74 * t66; t71 * t72, -t71 * t70, -t73, t71 * t65 + t73 * t69 + pkin(5) + qJ(2) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end