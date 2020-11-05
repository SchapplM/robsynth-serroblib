% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S4PRRR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 19:39
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S4PRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:17
	% EndTime: 2020-11-04 19:39:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:17
	% EndTime: 2020-11-04 19:39:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t68 = cos(pkin(8));
	t67 = sin(pkin(8));
	t1 = [t68, -t67, 0, 0; t67, t68, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:17
	% EndTime: 2020-11-04 19:39:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t69 = sin(pkin(8));
	t70 = sin(pkin(4));
	t78 = t69 * t70;
	t71 = cos(pkin(8));
	t77 = t71 * t70;
	t72 = cos(pkin(4));
	t73 = sin(qJ(2));
	t76 = t72 * t73;
	t74 = cos(qJ(2));
	t75 = t72 * t74;
	t1 = [-t69 * t76 + t71 * t74, -t69 * t75 - t71 * t73, t78, t71 * pkin(1) + pkin(5) * t78 + 0; t69 * t74 + t71 * t76, -t69 * t73 + t71 * t75, -t77, t69 * pkin(1) - pkin(5) * t77 + 0; t70 * t73, t70 * t74, t72, t72 * pkin(5) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:17
	% EndTime: 2020-11-04 19:39:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t82 = sin(pkin(4));
	t95 = t82 * pkin(5);
	t81 = sin(pkin(8));
	t84 = cos(pkin(4));
	t94 = t81 * t84;
	t85 = sin(qJ(3));
	t93 = t82 * t85;
	t87 = cos(qJ(3));
	t92 = t82 * t87;
	t83 = cos(pkin(8));
	t91 = t83 * t84;
	t86 = sin(qJ(2));
	t90 = t84 * t86;
	t88 = cos(qJ(2));
	t89 = t84 * t88;
	t80 = -t81 * t90 + t83 * t88;
	t79 = t81 * t88 + t83 * t90;
	t1 = [t80 * t87 + t81 * t93, -t80 * t85 + t81 * t92, t81 * t89 + t83 * t86, (t83 * pkin(2) + pkin(6) * t94) * t88 + (-pkin(2) * t94 + t83 * pkin(6)) * t86 + t81 * t95 + t83 * pkin(1) + 0; t79 * t87 - t83 * t93, -t79 * t85 - t83 * t92, t81 * t86 - t83 * t89, (t81 * pkin(2) - pkin(6) * t91) * t88 + (pkin(2) * t91 + t81 * pkin(6)) * t86 - t83 * t95 + t81 * pkin(1) + 0; t84 * t85 + t86 * t92, t84 * t87 - t86 * t93, -t82 * t88, t84 * pkin(5) + qJ(1) + 0 + (pkin(2) * t86 - pkin(6) * t88) * t82; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 19:39:17
	% EndTime: 2020-11-04 19:39:18
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (57->32), mult. (136->56), div. (0->0), fcn. (170->10), ass. (0->29)
	t103 = sin(pkin(4));
	t105 = cos(pkin(4));
	t108 = sin(qJ(2));
	t111 = cos(qJ(2));
	t107 = sin(qJ(3));
	t110 = cos(qJ(3));
	t115 = t110 * pkin(3) + t107 * pkin(7) + pkin(2);
	t113 = -pkin(6) * t111 + t115 * t108;
	t114 = t107 * pkin(3) - t110 * pkin(7) + pkin(5);
	t126 = t114 * t103 - t113 * t105;
	t122 = t103 * t110;
	t121 = t103 * t111;
	t120 = t105 * t108;
	t119 = t105 * t111;
	t118 = t107 * t103;
	t117 = t108 * t110;
	t116 = t110 * t111;
	t112 = pkin(6) * t108 + t115 * t111 + pkin(1);
	t109 = cos(qJ(4));
	t106 = sin(qJ(4));
	t104 = cos(pkin(8));
	t102 = sin(pkin(8));
	t101 = t105 * t117 - t118;
	t100 = t103 * t117 + t105 * t107;
	t99 = t102 * t119 + t104 * t108;
	t98 = t102 * t108 - t104 * t119;
	t97 = -t102 * t101 + t104 * t116;
	t96 = t104 * t101 + t102 * t116;
	t1 = [t99 * t106 + t97 * t109, -t97 * t106 + t99 * t109, -(t102 * t120 - t104 * t111) * t107 - t102 * t122, t126 * t102 + t112 * t104 + 0; t98 * t106 + t96 * t109, -t96 * t106 + t98 * t109, (t102 * t111 + t104 * t120) * t107 + t104 * t122, t112 * t102 - t126 * t104 + 0; t100 * t109 - t106 * t121, -t100 * t106 - t109 * t121, -t105 * t110 + t108 * t118, t113 * t103 + t114 * t105 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end