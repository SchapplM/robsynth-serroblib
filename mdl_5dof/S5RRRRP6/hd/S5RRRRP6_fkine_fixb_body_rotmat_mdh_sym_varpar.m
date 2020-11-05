% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRRRP6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:46
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:08
	% EndTime: 2020-11-04 20:46:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:08
	% EndTime: 2020-11-04 20:46:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t46 = cos(qJ(1));
	t45 = sin(qJ(1));
	t1 = [t46, -t45, 0, 0; t45, t46, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:08
	% EndTime: 2020-11-04 20:46:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t50 = cos(qJ(1));
	t49 = cos(qJ(2));
	t48 = sin(qJ(1));
	t47 = sin(qJ(2));
	t1 = [t50 * t49, -t50 * t47, t48, t50 * pkin(1) + t48 * pkin(6) + 0; t48 * t49, -t48 * t47, -t50, t48 * pkin(1) - t50 * pkin(6) + 0; t47, t49, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:08
	% EndTime: 2020-11-04 20:46:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t57 = pkin(7) + pkin(6);
	t56 = cos(qJ(1));
	t55 = sin(qJ(1));
	t54 = qJ(2) + qJ(3);
	t53 = cos(t54);
	t52 = sin(t54);
	t51 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t56 * t53, -t56 * t52, t55, t56 * t51 + t57 * t55 + 0; t55 * t53, -t55 * t52, -t56, t55 * t51 - t56 * t57 + 0; t52, t53, 0, sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:08
	% EndTime: 2020-11-04 20:46:08
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t62 = sin(qJ(4));
	t63 = sin(qJ(1));
	t71 = t63 * t62;
	t64 = cos(qJ(4));
	t70 = t63 * t64;
	t65 = cos(qJ(1));
	t69 = t65 * t62;
	t68 = t65 * t64;
	t61 = qJ(2) + qJ(3);
	t59 = sin(t61);
	t60 = cos(t61);
	t67 = pkin(3) * t60 + pkin(8) * t59 + cos(qJ(2)) * pkin(2) + pkin(1);
	t66 = pkin(7) + pkin(6);
	t1 = [t60 * t68 + t71, -t60 * t69 + t70, t65 * t59, t66 * t63 + t65 * t67 + 0; t60 * t70 - t69, -t60 * t71 - t68, t63 * t59, t63 * t67 - t65 * t66 + 0; t59 * t64, -t59 * t62, -t60, t59 * pkin(3) - t60 * pkin(8) + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:46:08
	% EndTime: 2020-11-04 20:46:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (45->22), mult. (40->24), div. (0->0), fcn. (53->8), ass. (0->16)
	t78 = sin(qJ(4));
	t79 = sin(qJ(1));
	t88 = t79 * t78;
	t80 = cos(qJ(4));
	t87 = t79 * t80;
	t81 = cos(qJ(1));
	t86 = t81 * t78;
	t85 = t81 * t80;
	t84 = pkin(4) * t78 + pkin(6) + pkin(7);
	t72 = t80 * pkin(4) + pkin(3);
	t76 = qJ(2) + qJ(3);
	t74 = sin(t76);
	t75 = cos(t76);
	t77 = -qJ(5) - pkin(8);
	t83 = t72 * t75 - t74 * t77 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t75 * t85 + t88, -t75 * t86 + t87, t81 * t74, t84 * t79 + t83 * t81 + 0; t75 * t87 - t86, -t75 * t88 - t85, t79 * t74, t83 * t79 - t84 * t81 + 0; t74 * t80, -t74 * t78, -t75, t74 * t72 + t75 * t77 + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end