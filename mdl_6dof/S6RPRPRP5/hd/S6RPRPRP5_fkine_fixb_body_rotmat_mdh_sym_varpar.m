% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRP5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:37
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:09
	% EndTime: 2020-11-04 21:37:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:09
	% EndTime: 2020-11-04 21:37:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t61 = cos(qJ(1));
	t60 = sin(qJ(1));
	t1 = [t61, -t60, 0, 0; t60, t61, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:09
	% EndTime: 2020-11-04 21:37:09
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t65 = cos(qJ(1));
	t64 = sin(qJ(1));
	t63 = cos(pkin(9));
	t62 = sin(pkin(9));
	t1 = [t65 * t63, -t65 * t62, t64, t65 * pkin(1) + t64 * qJ(2) + 0; t64 * t63, -t64 * t62, -t65, t64 * pkin(1) - t65 * qJ(2) + 0; t62, t63, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:09
	% EndTime: 2020-11-04 21:37:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t72 = cos(qJ(1));
	t71 = sin(qJ(1));
	t70 = pkin(7) + qJ(2);
	t69 = pkin(9) + qJ(3);
	t68 = cos(t69);
	t67 = sin(t69);
	t66 = cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t72 * t68, -t72 * t67, t71, t72 * t66 + t70 * t71 + 0; t71 * t68, -t71 * t67, -t72, t71 * t66 - t72 * t70 + 0; t67, t68, 0, sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:09
	% EndTime: 2020-11-04 21:37:09
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (37->19), mult. (33->22), div. (0->0), fcn. (46->8), ass. (0->14)
	t77 = sin(pkin(10));
	t80 = sin(qJ(1));
	t86 = t80 * t77;
	t78 = cos(pkin(10));
	t85 = t80 * t78;
	t81 = cos(qJ(1));
	t84 = t81 * t77;
	t83 = t81 * t78;
	t76 = pkin(9) + qJ(3);
	t74 = sin(t76);
	t75 = cos(t76);
	t82 = pkin(3) * t75 + qJ(4) * t74 + cos(pkin(9)) * pkin(2) + pkin(1);
	t79 = pkin(7) + qJ(2);
	t1 = [t75 * t83 + t86, -t75 * t84 + t85, t81 * t74, t79 * t80 + t82 * t81 + 0; t75 * t85 - t84, -t75 * t86 - t83, t80 * t74, -t81 * t79 + t82 * t80 + 0; t74 * t78, -t74 * t77, -t75, t74 * pkin(3) - t75 * qJ(4) + sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:09
	% EndTime: 2020-11-04 21:37:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (55->23), mult. (40->24), div. (0->0), fcn. (53->10), ass. (0->17)
	t93 = pkin(10) + qJ(5);
	t89 = sin(t93);
	t98 = sin(qJ(1));
	t105 = t98 * t89;
	t91 = cos(t93);
	t104 = t98 * t91;
	t99 = cos(qJ(1));
	t103 = t99 * t89;
	t102 = t99 * t91;
	t101 = sin(pkin(10)) * pkin(4) + pkin(7) + qJ(2);
	t87 = cos(pkin(10)) * pkin(4) + pkin(3);
	t94 = pkin(9) + qJ(3);
	t90 = sin(t94);
	t92 = cos(t94);
	t96 = -pkin(8) - qJ(4);
	t100 = t87 * t92 - t90 * t96 + cos(pkin(9)) * pkin(2) + pkin(1);
	t1 = [t92 * t102 + t105, -t92 * t103 + t104, t99 * t90, t100 * t99 + t101 * t98 + 0; t92 * t104 - t103, -t92 * t105 - t102, t98 * t90, t100 * t98 - t101 * t99 + 0; t90 * t91, -t90 * t89, -t92, t90 * t87 + t92 * t96 + sin(pkin(9)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:37:09
	% EndTime: 2020-11-04 21:37:09
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (80->28), mult. (60->30), div. (0->0), fcn. (77->10), ass. (0->21)
	t116 = pkin(10) + qJ(5);
	t112 = sin(t116);
	t121 = sin(qJ(1));
	t128 = t121 * t112;
	t114 = cos(t116);
	t127 = t121 * t114;
	t122 = cos(qJ(1));
	t126 = t122 * t112;
	t125 = t122 * t114;
	t124 = sin(pkin(10)) * pkin(4) + pkin(7) + qJ(2);
	t110 = cos(pkin(10)) * pkin(4) + pkin(3);
	t117 = pkin(9) + qJ(3);
	t113 = sin(t117);
	t115 = cos(t117);
	t119 = -pkin(8) - qJ(4);
	t123 = t110 * t115 - t113 * t119 + cos(pkin(9)) * pkin(2) + pkin(1);
	t109 = t115 * t125 + t128;
	t108 = t115 * t126 - t127;
	t107 = t115 * t127 - t126;
	t106 = t115 * t128 + t125;
	t1 = [t109, t122 * t113, t108, t109 * pkin(5) + t108 * qJ(6) + t124 * t121 + t123 * t122 + 0; t107, t121 * t113, t106, t107 * pkin(5) + t106 * qJ(6) + t123 * t121 - t124 * t122 + 0; t113 * t114, -t115, t113 * t112, sin(pkin(9)) * pkin(2) + t115 * t119 + pkin(6) + 0 + (pkin(5) * t114 + qJ(6) * t112 + t110) * t113; 0, 0, 0, 1;];
	Tc_mdh = t1;
end