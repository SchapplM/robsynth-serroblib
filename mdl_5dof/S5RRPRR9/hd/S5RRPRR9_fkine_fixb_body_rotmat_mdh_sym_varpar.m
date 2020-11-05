% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RRPRR9 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:37
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR9_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:39
	% EndTime: 2020-11-04 20:37:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:39
	% EndTime: 2020-11-04 20:37:39
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t53 = cos(qJ(1));
	t52 = sin(qJ(1));
	t1 = [t53, -t52, 0, 0; t52, t53, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:39
	% EndTime: 2020-11-04 20:37:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t57 = cos(qJ(1));
	t56 = cos(qJ(2));
	t55 = sin(qJ(1));
	t54 = sin(qJ(2));
	t1 = [t57 * t56, -t57 * t54, t55, t57 * pkin(1) + t55 * pkin(6) + 0; t55 * t56, -t55 * t54, -t57, t55 * pkin(1) - t57 * pkin(6) + 0; t54, t56, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:39
	% EndTime: 2020-11-04 20:37:39
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t64 = cos(qJ(1));
	t63 = sin(qJ(1));
	t62 = -qJ(3) - pkin(6);
	t61 = qJ(2) + pkin(9);
	t60 = cos(t61);
	t59 = sin(t61);
	t58 = cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t64 * t60, -t64 * t59, t63, t64 * t58 - t62 * t63 + 0; t63 * t60, -t63 * t59, -t64, t63 * t58 + t64 * t62 + 0; t59, t60, 0, sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:39
	% EndTime: 2020-11-04 20:37:39
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (37->21), mult. (35->25), div. (0->0), fcn. (48->10), ass. (0->17)
	t72 = sin(qJ(4));
	t74 = sin(qJ(1));
	t80 = t74 * t72;
	t75 = cos(qJ(4));
	t79 = t74 * t75;
	t76 = cos(qJ(1));
	t78 = t76 * t72;
	t77 = t76 * t75;
	t73 = sin(qJ(2));
	t71 = -qJ(3) - pkin(6);
	t70 = cos(pkin(9));
	t69 = sin(pkin(9));
	t68 = qJ(2) + pkin(9);
	t67 = cos(t68);
	t66 = sin(t68);
	t65 = (pkin(3) * t70 + pkin(7) * t69 + pkin(2)) * cos(qJ(2)) + (-pkin(3) * t69 + pkin(7) * t70) * t73 + pkin(1);
	t1 = [t67 * t77 + t80, -t67 * t78 + t79, t76 * t66, t65 * t76 - t71 * t74 + 0; t67 * t79 - t78, -t67 * t80 - t77, t74 * t66, t65 * t74 + t71 * t76 + 0; t66 * t75, -t66 * t72, -t67, pkin(2) * t73 + pkin(3) * t66 - pkin(7) * t67 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:37:39
	% EndTime: 2020-11-04 20:37:39
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t88 = qJ(4) + qJ(5);
	t85 = sin(t88);
	t91 = sin(qJ(1));
	t99 = t91 * t85;
	t86 = cos(t88);
	t98 = t91 * t86;
	t92 = cos(qJ(1));
	t97 = t92 * t85;
	t96 = t92 * t86;
	t95 = pkin(4) * sin(qJ(4)) + pkin(6) + qJ(3);
	t81 = cos(qJ(4)) * pkin(4) + pkin(3);
	t87 = qJ(2) + pkin(9);
	t83 = sin(t87);
	t84 = cos(t87);
	t93 = -pkin(8) - pkin(7);
	t94 = t81 * t84 - t83 * t93 + cos(qJ(2)) * pkin(2) + pkin(1);
	t1 = [t84 * t96 + t99, -t84 * t97 + t98, t92 * t83, t95 * t91 + t94 * t92 + 0; t84 * t98 - t97, -t84 * t99 - t96, t91 * t83, t94 * t91 - t95 * t92 + 0; t83 * t86, -t83 * t85, -t84, t83 * t81 + t84 * t93 + sin(qJ(2)) * pkin(2) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end