% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for the body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5RPPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 08:58:56
	% EndTime: 2022-01-23 08:58:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 08:58:56
	% EndTime: 2022-01-23 08:58:56
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t62 = cos(qJ(1));
	t61 = sin(qJ(1));
	t1 = [t62, -t61, 0, 0; t61, t62, 0, 0; 0, 0, 1, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 08:58:56
	% EndTime: 2022-01-23 08:58:56
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t66 = cos(qJ(1));
	t65 = sin(qJ(1));
	t64 = cos(pkin(7));
	t63 = sin(pkin(7));
	t1 = [t66 * t64, -t66 * t63, t65, t66 * pkin(1) + t65 * qJ(2) + 0; t65 * t64, -t65 * t63, -t66, t65 * pkin(1) - t66 * qJ(2) + 0; t63, t64, 0, pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 08:58:56
	% EndTime: 2022-01-23 08:58:56
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (17->15), mult. (26->20), div. (0->0), fcn. (39->6), ass. (0->12)
	t68 = sin(pkin(8));
	t72 = sin(qJ(1));
	t77 = t72 * t68;
	t70 = cos(pkin(8));
	t76 = t72 * t70;
	t73 = cos(qJ(1));
	t75 = t73 * t68;
	t74 = t73 * t70;
	t71 = cos(pkin(7));
	t69 = sin(pkin(7));
	t67 = pkin(2) * t71 + t69 * qJ(3) + pkin(1);
	t1 = [t71 * t74 + t77, -t71 * t75 + t76, t73 * t69, t72 * qJ(2) + t67 * t73 + 0; t71 * t76 - t75, -t71 * t77 - t74, t72 * t69, -t73 * qJ(2) + t67 * t72 + 0; t69 * t70, -t69 * t68, -t71, t69 * pkin(2) - t71 * qJ(3) + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 08:58:56
	% EndTime: 2022-01-23 08:58:56
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (33->23), mult. (57->36), div. (0->0), fcn. (78->8), ass. (0->21)
	t85 = sin(pkin(7));
	t87 = cos(pkin(8));
	t97 = t85 * t87;
	t89 = sin(qJ(1));
	t96 = t85 * t89;
	t90 = cos(qJ(1));
	t95 = t85 * t90;
	t84 = sin(pkin(8));
	t94 = t89 * t84;
	t93 = t89 * t87;
	t92 = t90 * t84;
	t91 = t90 * t87;
	t88 = cos(pkin(7));
	t86 = cos(pkin(9));
	t83 = sin(pkin(9));
	t82 = t87 * pkin(3) + t84 * qJ(4) + pkin(2);
	t81 = -t84 * pkin(3) + qJ(4) * t87 - qJ(2);
	t80 = t88 * t91 + t94;
	t79 = t88 * t93 - t92;
	t78 = t85 * qJ(3) + t82 * t88 + pkin(1);
	t1 = [t80 * t86 + t83 * t95, -t80 * t83 + t86 * t95, t88 * t92 - t93, t78 * t90 - t81 * t89 + 0; t79 * t86 + t83 * t96, -t79 * t83 + t86 * t96, t88 * t94 + t91, t78 * t89 + t81 * t90 + 0; -t88 * t83 + t86 * t97, -t83 * t97 - t88 * t86, t85 * t84, -t88 * qJ(3) + t82 * t85 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 08:58:56
	% EndTime: 2022-01-23 08:58:56
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (62->36), mult. (109->59), div. (0->0), fcn. (143->10), ass. (0->28)
	t106 = sin(pkin(8));
	t111 = sin(qJ(5));
	t124 = t106 * t111;
	t113 = cos(qJ(5));
	t123 = t106 * t113;
	t107 = sin(pkin(7));
	t108 = cos(pkin(9));
	t122 = t107 * t108;
	t109 = cos(pkin(8));
	t121 = t107 * t109;
	t110 = cos(pkin(7));
	t120 = t110 * t108;
	t112 = sin(qJ(1));
	t119 = t112 * t106;
	t114 = cos(qJ(1));
	t118 = t114 * t106;
	t117 = t114 * t109;
	t105 = sin(pkin(9));
	t116 = t105 * pkin(4) - t108 * pkin(6) + qJ(3);
	t104 = t108 * pkin(4) + t105 * pkin(6) + pkin(3);
	t115 = qJ(4) * t109 - t104 * t106 - qJ(2);
	t103 = t108 * t124 + t113 * t109;
	t102 = t107 * t105 + t109 * t120;
	t101 = -t110 * t105 + t108 * t121;
	t100 = t106 * qJ(4) + t104 * t109 + pkin(2);
	t99 = -t102 * t111 + t110 * t123;
	t98 = t100 * t110 + t116 * t107 + pkin(1);
	t1 = [(t102 * t113 + t110 * t124) * t114 + t112 * (t108 * t123 - t111 * t109), -t112 * t103 + t99 * t114, (t110 * t117 + t119) * t105 - t114 * t122, -t115 * t112 + t98 * t114 + 0; (t102 * t112 - t108 * t118) * t113 + (t110 * t119 + t117) * t111, t114 * t103 + t99 * t112, (t112 * t110 * t109 - t118) * t105 - t112 * t122, t98 * t112 + t115 * t114 + 0; t101 * t113 + t107 * t124, -t101 * t111 + t107 * t123, t105 * t121 + t120, t100 * t107 - t116 * t110 + pkin(5) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end