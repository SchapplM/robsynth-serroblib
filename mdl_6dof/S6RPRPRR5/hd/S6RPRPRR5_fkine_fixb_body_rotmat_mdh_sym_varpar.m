% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RPRPRR5 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:40
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RPRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:31
	% EndTime: 2020-11-04 21:40:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:31
	% EndTime: 2020-11-04 21:40:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t71 = cos(qJ(1));
	t70 = sin(qJ(1));
	t1 = [t71, -t70, 0, 0; t70, t71, 0, 0; 0, 0, 1, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:31
	% EndTime: 2020-11-04 21:40:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (8->8), mult. (8->8), div. (0->0), fcn. (16->4), ass. (0->5)
	t75 = cos(qJ(1));
	t74 = sin(qJ(1));
	t73 = cos(pkin(10));
	t72 = sin(pkin(10));
	t1 = [t75 * t73, -t75 * t72, t74, t75 * pkin(1) + t74 * qJ(2) + 0; t74 * t73, -t74 * t72, -t75, t74 * pkin(1) - t75 * qJ(2) + 0; t72, t73, 0, pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:31
	% EndTime: 2020-11-04 21:40:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (19->12), mult. (11->10), div. (0->0), fcn. (19->6), ass. (0->8)
	t82 = cos(qJ(1));
	t81 = sin(qJ(1));
	t80 = pkin(7) + qJ(2);
	t79 = pkin(10) + qJ(3);
	t78 = cos(t79);
	t77 = sin(t79);
	t76 = cos(pkin(10)) * pkin(2) + pkin(1);
	t1 = [t82 * t78, -t82 * t77, t81, t82 * t76 + t80 * t81 + 0; t81 * t78, -t81 * t77, -t82, t81 * t76 - t82 * t80 + 0; t77, t78, 0, sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:31
	% EndTime: 2020-11-04 21:40:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->15), mult. (21->14), div. (0->0), fcn. (29->6), ass. (0->8)
	t86 = pkin(10) + qJ(3);
	t84 = sin(t86);
	t85 = cos(t86);
	t90 = pkin(3) * t85 + qJ(4) * t84 + cos(pkin(10)) * pkin(2) + pkin(1);
	t89 = cos(qJ(1));
	t88 = sin(qJ(1));
	t87 = pkin(7) + qJ(2);
	t1 = [t89 * t85, t88, t89 * t84, t87 * t88 + t90 * t89 + 0; t88 * t85, -t89, t88 * t84, -t89 * t87 + t90 * t88 + 0; t84, 0, -t85, t84 * pkin(3) - t85 * qJ(4) + sin(pkin(10)) * pkin(2) + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:31
	% EndTime: 2020-11-04 21:40:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (48->21), mult. (37->22), div. (0->0), fcn. (51->10), ass. (0->15)
	t106 = cos(qJ(5));
	t105 = sin(qJ(5));
	t104 = pkin(3) + pkin(4);
	t103 = cos(qJ(1));
	t102 = sin(qJ(1));
	t101 = cos(pkin(10));
	t100 = sin(pkin(10));
	t99 = pkin(10) + qJ(3);
	t98 = qJ(2) + pkin(7) - pkin(8);
	t97 = cos(t99);
	t96 = sin(t99);
	t93 = -t97 * t105 + t96 * t106;
	t92 = -t96 * t105 - t97 * t106;
	t91 = (qJ(4) * t100 + t104 * t101) * cos(qJ(3)) + (qJ(4) * t101 - t100 * t104) * sin(qJ(3)) + t101 * pkin(2) + pkin(1);
	t1 = [-t103 * t92, t103 * t93, -t102, t98 * t102 + t91 * t103 + 0; -t102 * t92, t102 * t93, t103, t91 * t102 - t98 * t103 + 0; t93, t92, 0, t100 * pkin(2) - t97 * qJ(4) + t104 * t96 + pkin(6) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:40:31
	% EndTime: 2020-11-04 21:40:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (85->34), mult. (81->40), div. (0->0), fcn. (103->14), ass. (0->24)
	t119 = sin(qJ(6));
	t121 = sin(qJ(1));
	t129 = t121 * t119;
	t122 = cos(qJ(6));
	t128 = t121 * t122;
	t124 = cos(qJ(1));
	t127 = t124 * t119;
	t126 = t124 * t122;
	t116 = pkin(10) + qJ(3);
	t125 = pkin(3) + pkin(4);
	t123 = cos(qJ(5));
	t120 = sin(qJ(5));
	t118 = cos(pkin(10));
	t117 = sin(pkin(10));
	t115 = qJ(2) + pkin(7) - pkin(8);
	t114 = -qJ(5) + t116;
	t113 = cos(t116);
	t112 = sin(t116);
	t111 = t118 * pkin(5) - t117 * pkin(9);
	t110 = t117 * pkin(5) + t118 * pkin(9);
	t109 = t112 * t120 + t113 * t123;
	t108 = t112 * t123 - t113 * t120;
	t107 = (qJ(4) * t117 + t110 * t120 + t111 * t123 + t125 * t118) * cos(qJ(3)) + (qJ(4) * t118 - t110 * t123 + t111 * t120 - t117 * t125) * sin(qJ(3)) + t118 * pkin(2) + pkin(1);
	t1 = [t109 * t126 - t129, -t109 * t127 - t128, -t124 * t108, t107 * t124 + t115 * t121 + 0; t109 * t128 + t127, -t109 * t129 + t126, -t121 * t108, t107 * t121 - t115 * t124 + 0; t108 * t122, -t108 * t119, t109, pkin(9) * cos(t114) + pkin(5) * sin(t114) + t125 * t112 - t113 * qJ(4) + t117 * pkin(2) + 0 + pkin(6); 0, 0, 0, 1;];
	Tc_mdh = t1;
end