% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR6 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:00
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:50
	% EndTime: 2020-11-04 20:00:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:50
	% EndTime: 2020-11-04 20:00:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t73 = cos(pkin(9));
	t72 = sin(pkin(9));
	t1 = [t73, -t72, 0, 0; t72, t73, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:50
	% EndTime: 2020-11-04 20:00:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t74 = sin(pkin(9));
	t75 = sin(pkin(5));
	t83 = t74 * t75;
	t76 = cos(pkin(9));
	t82 = t76 * t75;
	t77 = cos(pkin(5));
	t78 = sin(qJ(2));
	t81 = t77 * t78;
	t79 = cos(qJ(2));
	t80 = t77 * t79;
	t1 = [-t74 * t81 + t76 * t79, -t74 * t80 - t76 * t78, t83, t76 * pkin(1) + pkin(6) * t83 + 0; t74 * t79 + t76 * t81, -t74 * t78 + t76 * t80, -t82, t74 * pkin(1) - pkin(6) * t82 + 0; t75 * t78, t75 * t79, t77, t77 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:50
	% EndTime: 2020-11-04 20:00:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->46), div. (0->0), fcn. (85->8), ass. (0->20)
	t87 = sin(pkin(9));
	t102 = t87 * pkin(2);
	t90 = cos(pkin(9));
	t101 = t90 * pkin(2);
	t88 = sin(pkin(5));
	t100 = t87 * t88;
	t99 = t88 * t90;
	t92 = sin(qJ(2));
	t98 = t88 * t92;
	t91 = cos(pkin(5));
	t97 = t91 * t92;
	t93 = cos(qJ(2));
	t96 = t91 * t93;
	t95 = t87 * qJ(3);
	t94 = t90 * qJ(3);
	t89 = cos(pkin(10));
	t86 = sin(pkin(10));
	t85 = t87 * t93 + t90 * t97;
	t84 = t87 * t97 - t90 * t93;
	t1 = [t86 * t100 - t84 * t89, t89 * t100 + t84 * t86, t87 * t96 + t90 * t92, (t91 * t95 + t101) * t93 + (-t91 * t102 + t94) * t92 + pkin(6) * t100 + t90 * pkin(1) + 0; t85 * t89 - t86 * t99, -t85 * t86 - t89 * t99, t87 * t92 - t90 * t96, (-t91 * t94 + t102) * t93 + (t91 * t101 + t95) * t92 - pkin(6) * t99 + t87 * pkin(1) + 0; t91 * t86 + t89 * t98, -t86 * t98 + t91 * t89, -t88 * t93, t91 * pkin(6) + qJ(1) + 0 + (pkin(2) * t92 - qJ(3) * t93) * t88; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:50
	% EndTime: 2020-11-04 20:00:50
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (54->31), mult. (72->49), div. (0->0), fcn. (93->10), ass. (0->23)
	t106 = cos(pkin(10)) * pkin(3) + pkin(2);
	t110 = sin(pkin(9));
	t124 = t110 * t106;
	t111 = sin(pkin(5));
	t123 = t110 * t111;
	t112 = cos(pkin(9));
	t122 = t111 * t112;
	t115 = sin(qJ(2));
	t121 = t111 * t115;
	t120 = t112 * t106;
	t113 = cos(pkin(5));
	t114 = qJ(3) + pkin(7);
	t119 = t113 * t114;
	t118 = t113 * t115;
	t116 = cos(qJ(2));
	t117 = t113 * t116;
	t109 = pkin(10) + qJ(4);
	t108 = cos(t109);
	t107 = sin(t109);
	t105 = sin(pkin(10)) * pkin(3) + pkin(6);
	t104 = t110 * t116 + t112 * t118;
	t103 = t110 * t118 - t112 * t116;
	t1 = [-t103 * t108 + t107 * t123, t103 * t107 + t108 * t123, t110 * t117 + t112 * t115, (t110 * t119 + t120) * t116 + (t112 * t114 - t113 * t124) * t115 + t105 * t123 + t112 * pkin(1) + 0; t104 * t108 - t107 * t122, -t104 * t107 - t108 * t122, t110 * t115 - t112 * t117, (-t112 * t119 + t124) * t116 + (t110 * t114 + t113 * t120) * t115 - t105 * t122 + t110 * pkin(1) + 0; t113 * t107 + t108 * t121, -t107 * t121 + t113 * t108, -t111 * t116, t105 * t113 + qJ(1) + 0 + (t106 * t115 - t114 * t116) * t111; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:00:50
	% EndTime: 2020-11-04 20:00:50
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (100->43), mult. (142->67), div. (0->0), fcn. (186->12), ass. (0->34)
	t136 = cos(pkin(10)) * pkin(3) + pkin(2);
	t140 = sin(pkin(9));
	t157 = t140 * t136;
	t141 = sin(pkin(5));
	t156 = t140 * t141;
	t142 = cos(pkin(9));
	t155 = t141 * t142;
	t146 = sin(qJ(2));
	t154 = t141 * t146;
	t148 = cos(qJ(2));
	t153 = t141 * t148;
	t152 = t142 * t136;
	t143 = cos(pkin(5));
	t144 = qJ(3) + pkin(7);
	t151 = t143 * t144;
	t150 = t143 * t146;
	t149 = t143 * t148;
	t147 = cos(qJ(5));
	t145 = sin(qJ(5));
	t139 = pkin(10) + qJ(4);
	t138 = cos(t139);
	t137 = sin(t139);
	t135 = sin(pkin(10)) * pkin(3) + pkin(6);
	t134 = t140 * t149 + t142 * t146;
	t133 = t140 * t148 + t142 * t150;
	t132 = t140 * t146 - t142 * t149;
	t131 = t140 * t150 - t142 * t148;
	t130 = t143 * t137 + t138 * t154;
	t129 = t137 * t154 - t143 * t138;
	t128 = -t131 * t138 + t137 * t156;
	t127 = t133 * t138 - t137 * t155;
	t126 = t133 * t137 + t138 * t155;
	t125 = t131 * t137 + t138 * t156;
	t1 = [t128 * t147 + t134 * t145, -t128 * t145 + t134 * t147, -t125, t128 * pkin(4) - t125 * pkin(8) + (t140 * t151 + t152) * t148 + (t142 * t144 - t143 * t157) * t146 + t135 * t156 + t142 * pkin(1) + 0; t127 * t147 + t132 * t145, -t127 * t145 + t132 * t147, t126, t127 * pkin(4) + t126 * pkin(8) + (-t142 * t151 + t157) * t148 + (t140 * t144 + t143 * t152) * t146 - t135 * t155 + t140 * pkin(1) + 0; t130 * t147 - t145 * t153, -t130 * t145 - t147 * t153, t129, t130 * pkin(4) + t129 * pkin(8) + t135 * t143 + qJ(1) + 0 + (t136 * t146 - t144 * t148) * t141; 0, 0, 0, 1;];
	Tc_mdh = t1;
end