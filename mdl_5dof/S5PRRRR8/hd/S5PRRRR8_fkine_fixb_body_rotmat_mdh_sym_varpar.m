% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR8 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 20:09
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S5PRRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR8_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:23
	% EndTime: 2020-11-04 20:09:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:23
	% EndTime: 2020-11-04 20:09:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t71 = cos(pkin(10));
	t70 = sin(pkin(10));
	t1 = [t71, -t70, 0, 0; t70, t71, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:23
	% EndTime: 2020-11-04 20:09:23
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t72 = sin(pkin(10));
	t73 = sin(pkin(5));
	t81 = t72 * t73;
	t74 = cos(pkin(10));
	t80 = t74 * t73;
	t75 = cos(pkin(5));
	t76 = sin(qJ(2));
	t79 = t75 * t76;
	t77 = cos(qJ(2));
	t78 = t75 * t77;
	t1 = [-t72 * t79 + t74 * t77, -t72 * t78 - t74 * t76, t81, t74 * pkin(1) + pkin(6) * t81 + 0; t72 * t77 + t74 * t79, -t72 * t76 + t74 * t78, -t80, t72 * pkin(1) - pkin(6) * t80 + 0; t73 * t76, t73 * t77, t75, t75 * pkin(6) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:23
	% EndTime: 2020-11-04 20:09:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (29->27), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->18)
	t85 = sin(pkin(5));
	t98 = t85 * pkin(6);
	t84 = sin(pkin(10));
	t87 = cos(pkin(5));
	t97 = t84 * t87;
	t88 = sin(qJ(3));
	t96 = t85 * t88;
	t90 = cos(qJ(3));
	t95 = t85 * t90;
	t86 = cos(pkin(10));
	t94 = t86 * t87;
	t89 = sin(qJ(2));
	t93 = t87 * t89;
	t91 = cos(qJ(2));
	t92 = t87 * t91;
	t83 = t84 * t91 + t86 * t93;
	t82 = t84 * t93 - t86 * t91;
	t1 = [-t82 * t90 + t84 * t96, t82 * t88 + t84 * t95, t84 * t92 + t86 * t89, (t86 * pkin(2) + pkin(7) * t97) * t91 + (-pkin(2) * t97 + t86 * pkin(7)) * t89 + t84 * t98 + t86 * pkin(1) + 0; t83 * t90 - t86 * t96, -t83 * t88 - t86 * t95, t84 * t89 - t86 * t92, (t84 * pkin(2) - pkin(7) * t94) * t91 + (pkin(2) * t94 + t84 * pkin(7)) * t89 - t86 * t98 + t84 * pkin(1) + 0; t87 * t88 + t89 * t95, t87 * t90 - t89 * t96, -t85 * t91, t87 * pkin(6) + qJ(1) + 0 + (pkin(2) * t89 - pkin(7) * t91) * t85; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:23
	% EndTime: 2020-11-04 20:09:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (50->27), mult. (75->40), div. (0->0), fcn. (100->10), ass. (0->22)
	t108 = sin(pkin(10));
	t109 = sin(pkin(5));
	t121 = t108 * t109;
	t110 = cos(pkin(10));
	t120 = t109 * t110;
	t113 = sin(qJ(2));
	t119 = t109 * t113;
	t111 = cos(pkin(5));
	t118 = t111 * t113;
	t114 = cos(qJ(2));
	t117 = t111 * t114;
	t116 = pkin(3) * sin(qJ(3)) + pkin(6);
	t115 = pkin(7) + pkin(8);
	t107 = qJ(3) + qJ(4);
	t106 = cos(t107);
	t105 = sin(t107);
	t104 = cos(qJ(3)) * pkin(3) + pkin(2);
	t102 = t108 * t117 + t110 * t113;
	t101 = t108 * t114 + t110 * t118;
	t100 = t108 * t113 - t110 * t117;
	t99 = t108 * t118 - t110 * t114;
	t1 = [t105 * t121 - t99 * t106, t99 * t105 + t106 * t121, t102, t110 * pkin(1) + t102 * t115 - t99 * t104 + t116 * t121 + 0; t101 * t106 - t105 * t120, -t101 * t105 - t106 * t120, t100, t108 * pkin(1) + t100 * t115 + t101 * t104 - t116 * t120 + 0; t111 * t105 + t106 * t119, -t105 * t119 + t111 * t106, -t109 * t114, qJ(1) + 0 + t116 * t111 + (t104 * t113 - t114 * t115) * t109; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 20:09:23
	% EndTime: 2020-11-04 20:09:23
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (96->39), mult. (145->58), div. (0->0), fcn. (193->12), ass. (0->31)
	t137 = sin(pkin(10));
	t138 = sin(pkin(5));
	t153 = t137 * t138;
	t139 = cos(pkin(10));
	t152 = t138 * t139;
	t143 = sin(qJ(2));
	t151 = t138 * t143;
	t145 = cos(qJ(2));
	t150 = t138 * t145;
	t140 = cos(pkin(5));
	t149 = t140 * t143;
	t148 = t140 * t145;
	t147 = pkin(3) * sin(qJ(3)) + pkin(6);
	t146 = pkin(7) + pkin(8);
	t144 = cos(qJ(5));
	t141 = sin(qJ(5));
	t136 = qJ(3) + qJ(4);
	t135 = cos(t136);
	t134 = sin(t136);
	t133 = cos(qJ(3)) * pkin(3) + pkin(2);
	t131 = t137 * t148 + t139 * t143;
	t130 = t137 * t145 + t139 * t149;
	t129 = t137 * t143 - t139 * t148;
	t128 = t137 * t149 - t139 * t145;
	t127 = t140 * t134 + t135 * t151;
	t126 = t134 * t151 - t140 * t135;
	t125 = -t128 * t135 + t134 * t153;
	t124 = t130 * t135 - t134 * t152;
	t123 = t130 * t134 + t135 * t152;
	t122 = t128 * t134 + t135 * t153;
	t1 = [t125 * t144 + t131 * t141, -t125 * t141 + t131 * t144, -t122, t139 * pkin(1) + t125 * pkin(4) - t122 * pkin(9) - t128 * t133 + t131 * t146 + t147 * t153 + 0; t124 * t144 + t129 * t141, -t124 * t141 + t129 * t144, t123, t137 * pkin(1) + t124 * pkin(4) + t123 * pkin(9) + t129 * t146 + t130 * t133 - t147 * t152 + 0; t127 * t144 - t141 * t150, -t127 * t141 - t144 * t150, t126, t127 * pkin(4) + t126 * pkin(9) + qJ(1) + 0 + t147 * t140 + (t133 * t143 - t145 * t146) * t138; 0, 0, 0, 1;];
	Tc_mdh = t1;
end