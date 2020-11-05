% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PPRRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:55
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:33
	% EndTime: 2020-11-04 19:55:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:33
	% EndTime: 2020-11-04 19:55:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t54 = cos(pkin(8));
	t53 = sin(pkin(8));
	t1 = [t54, -t53, 0, 0; t53, t54, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:33
	% EndTime: 2020-11-04 19:55:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t58 = cos(pkin(8));
	t57 = cos(pkin(9));
	t56 = sin(pkin(8));
	t55 = sin(pkin(9));
	t1 = [t58 * t57, -t58 * t55, t56, t58 * pkin(1) + t56 * qJ(2) + 0; t56 * t57, -t56 * t55, -t58, t56 * pkin(1) - t58 * qJ(2) + 0; t55, t57, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:33
	% EndTime: 2020-11-04 19:55:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->12)
	t61 = sin(pkin(8));
	t64 = sin(qJ(3));
	t69 = t61 * t64;
	t65 = cos(qJ(3));
	t68 = t61 * t65;
	t63 = cos(pkin(8));
	t67 = t63 * t64;
	t66 = t63 * t65;
	t62 = cos(pkin(9));
	t60 = sin(pkin(9));
	t59 = t62 * pkin(2) + t60 * pkin(5) + pkin(1);
	t1 = [t62 * t66 + t69, -t62 * t67 + t68, t63 * t60, t61 * qJ(2) + t59 * t63 + 0; t62 * t68 - t67, -t62 * t69 - t66, t61 * t60, -t63 * qJ(2) + t59 * t61 + 0; t60 * t65, -t60 * t64, -t62, t60 * pkin(2) - t62 * pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:33
	% EndTime: 2020-11-04 19:55:33
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (33->29), mult. (67->48), div. (0->0), fcn. (88->8), ass. (0->21)
	t73 = sin(pkin(9));
	t77 = sin(qJ(4));
	t89 = t73 * t77;
	t79 = cos(qJ(4));
	t88 = t73 * t79;
	t80 = cos(qJ(3));
	t87 = t73 * t80;
	t74 = sin(pkin(8));
	t75 = cos(pkin(9));
	t86 = t74 * t75;
	t78 = sin(qJ(3));
	t85 = t74 * t78;
	t84 = t74 * t80;
	t76 = cos(pkin(8));
	t83 = t75 * t76;
	t82 = t76 * t78;
	t81 = t76 * t80;
	t72 = t75 * pkin(2) + t73 * pkin(5) + pkin(1);
	t71 = t75 * t81 + t85;
	t70 = t75 * t84 - t82;
	t1 = [t71 * t79 + t76 * t89, -t71 * t77 + t76 * t88, t75 * t82 - t84, (pkin(3) * t83 - t74 * pkin(6)) * t80 + (t74 * pkin(3) + pkin(6) * t83) * t78 + t72 * t76 + t74 * qJ(2) + 0; t70 * t79 + t74 * t89, -t70 * t77 + t74 * t88, t75 * t85 + t81, (pkin(3) * t86 + t76 * pkin(6)) * t80 + (-t76 * pkin(3) + pkin(6) * t86) * t78 + t72 * t74 - t76 * qJ(2) + 0; -t75 * t77 + t79 * t87, -t75 * t79 - t77 * t87, t73 * t78, -t75 * pkin(5) + qJ(1) + 0 + (pkin(3) * t80 + pkin(6) * t78 + pkin(2)) * t73; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:55:33
	% EndTime: 2020-11-04 19:55:33
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (54->29), mult. (78->40), div. (0->0), fcn. (103->10), ass. (0->25)
	t116 = pkin(5) + pkin(4) * sin(qJ(4));
	t100 = sin(pkin(8));
	t99 = sin(pkin(9));
	t114 = t100 * t99;
	t102 = cos(pkin(8));
	t113 = t102 * t99;
	t105 = cos(qJ(3));
	t112 = t105 * t99;
	t104 = sin(qJ(3));
	t111 = t100 * t104;
	t110 = t100 * t105;
	t109 = t102 * t104;
	t108 = t102 * t105;
	t101 = cos(pkin(9));
	t107 = t101 * pkin(2) + t116 * t99 + pkin(1);
	t106 = pkin(6) + pkin(7);
	t98 = qJ(4) + qJ(5);
	t97 = cos(t98);
	t96 = sin(t98);
	t95 = cos(qJ(4)) * pkin(4) + pkin(3);
	t93 = t101 * t108 + t111;
	t92 = t101 * t109 - t110;
	t91 = t101 * t110 - t109;
	t90 = t101 * t111 + t108;
	t1 = [t96 * t113 + t93 * t97, t97 * t113 - t93 * t96, t92, t100 * qJ(2) + t107 * t102 + t92 * t106 + t93 * t95 + 0; t96 * t114 + t91 * t97, t97 * t114 - t91 * t96, t90, -t102 * qJ(2) + t107 * t100 + t90 * t106 + t91 * t95 + 0; -t101 * t96 + t97 * t112, -t101 * t97 - t96 * t112, t99 * t104, qJ(1) + 0 - t116 * t101 + (t104 * t106 + t105 * t95 + pkin(2)) * t99; 0, 0, 0, 1;];
	Tc_mdh = t1;
end