% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRPP3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:02
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:02:30
	% EndTime: 2020-11-04 20:02:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:02:30
	% EndTime: 2020-11-04 20:02:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t59 = cos(pkin(7));
	t58 = sin(pkin(7));
	t1 = [t59, -t58, 0, 0; t58, t59, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:02:30
	% EndTime: 2020-11-04 20:02:30
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t63 = cos(qJ(2));
	t62 = sin(qJ(2));
	t61 = cos(pkin(7));
	t60 = sin(pkin(7));
	t1 = [t61 * t63, -t61 * t62, t60, t61 * pkin(1) + t60 * pkin(5) + 0; t60 * t63, -t60 * t62, -t61, t60 * pkin(1) - t61 * pkin(5) + 0; t62, t63, 0, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:02:30
	% EndTime: 2020-11-04 20:02:30
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (30->22), div. (0->0), fcn. (43->6), ass. (0->10)
	t66 = sin(qJ(3));
	t69 = cos(qJ(2));
	t72 = t66 * t69;
	t68 = cos(qJ(3));
	t71 = t68 * t69;
	t67 = sin(qJ(2));
	t70 = pkin(2) * t69 + pkin(6) * t67 + pkin(1);
	t65 = cos(pkin(7));
	t64 = sin(pkin(7));
	t1 = [t64 * t66 + t65 * t71, t64 * t68 - t65 * t72, t65 * t67, t64 * pkin(5) + t70 * t65 + 0; t64 * t71 - t65 * t66, -t64 * t72 - t65 * t68, t64 * t67, -t65 * pkin(5) + t70 * t64 + 0; t67 * t68, -t67 * t66, -t69, t67 * pkin(2) - t69 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:02:30
	% EndTime: 2020-11-04 20:02:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (33->23), mult. (65->38), div. (0->0), fcn. (86->8), ass. (0->19)
	t77 = sin(pkin(7));
	t81 = sin(qJ(2));
	t90 = t77 * t81;
	t79 = cos(pkin(7));
	t89 = t79 * t81;
	t80 = sin(qJ(3));
	t83 = cos(qJ(2));
	t88 = t80 * t83;
	t82 = cos(qJ(3));
	t87 = t81 * t82;
	t86 = t82 * t83;
	t75 = t82 * pkin(3) + t80 * qJ(4) + pkin(2);
	t85 = pkin(6) * t81 + t75 * t83 + pkin(1);
	t84 = pkin(3) * t80 - qJ(4) * t82 + pkin(5);
	t78 = cos(pkin(8));
	t76 = sin(pkin(8));
	t74 = t77 * t80 + t79 * t86;
	t73 = t77 * t86 - t79 * t80;
	t1 = [t74 * t78 + t76 * t89, -t74 * t76 + t78 * t89, -t77 * t82 + t79 * t88, t84 * t77 + t85 * t79 + 0; t73 * t78 + t76 * t90, -t73 * t76 + t78 * t90, t77 * t88 + t79 * t82, t85 * t77 - t84 * t79 + 0; -t83 * t76 + t78 * t87, -t76 * t87 - t83 * t78, t81 * t80, -t83 * pkin(6) + t75 * t81 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:02:31
	% EndTime: 2020-11-04 20:02:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (49->27), mult. (81->42), div. (0->0), fcn. (102->8), ass. (0->21)
	t101 = sin(qJ(2));
	t97 = sin(pkin(7));
	t110 = t101 * t97;
	t99 = cos(pkin(7));
	t109 = t101 * t99;
	t100 = sin(qJ(3));
	t103 = cos(qJ(2));
	t108 = t100 * t103;
	t102 = cos(qJ(3));
	t107 = t101 * t102;
	t106 = t102 * t103;
	t96 = sin(pkin(8));
	t98 = cos(pkin(8));
	t94 = t98 * pkin(4) + qJ(5) * t96 + pkin(3);
	t91 = t100 * qJ(4) + t94 * t102 + pkin(2);
	t95 = -t96 * pkin(4) + qJ(5) * t98 - pkin(6);
	t105 = -t101 * t95 + t103 * t91 + pkin(1);
	t104 = qJ(4) * t102 - t100 * t94 - pkin(5);
	t93 = t97 * t100 + t99 * t106;
	t92 = -t99 * t100 + t97 * t106;
	t1 = [t96 * t109 + t93 * t98, -t97 * t102 + t99 * t108, -t98 * t109 + t93 * t96, -t104 * t97 + t105 * t99 + 0; t96 * t110 + t92 * t98, t99 * t102 + t97 * t108, -t98 * t110 + t92 * t96, t104 * t99 + t105 * t97 + 0; -t103 * t96 + t98 * t107, t101 * t100, t103 * t98 + t96 * t107, t91 * t101 + t95 * t103 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end