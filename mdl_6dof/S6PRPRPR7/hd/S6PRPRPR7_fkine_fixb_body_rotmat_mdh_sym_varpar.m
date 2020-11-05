% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR7 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:01
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:02
	% EndTime: 2020-11-04 21:01:02
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:02
	% EndTime: 2020-11-04 21:01:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t98 = cos(pkin(10));
	t97 = sin(pkin(10));
	t1 = [t98, -t97, 0, 0; t97, t98, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:02
	% EndTime: 2020-11-04 21:01:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t100 = sin(pkin(6));
	t99 = sin(pkin(10));
	t108 = t99 * t100;
	t101 = cos(pkin(10));
	t107 = t101 * t100;
	t102 = cos(pkin(6));
	t103 = sin(qJ(2));
	t106 = t102 * t103;
	t104 = cos(qJ(2));
	t105 = t102 * t104;
	t1 = [t101 * t104 - t99 * t106, -t101 * t103 - t99 * t105, t108, t101 * pkin(1) + pkin(7) * t108 + 0; t101 * t106 + t99 * t104, t101 * t105 - t99 * t103, -t107, t99 * pkin(1) - pkin(7) * t107 + 0; t100 * t103, t100 * t104, t102, t102 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:02
	% EndTime: 2020-11-04 21:01:02
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (23->23), mult. (43->36), div. (0->0), fcn. (56->6), ass. (0->13)
	t109 = sin(pkin(10));
	t110 = sin(pkin(6));
	t120 = t109 * t110;
	t112 = cos(pkin(6));
	t119 = t109 * t112;
	t111 = cos(pkin(10));
	t118 = t111 * t110;
	t117 = t111 * t112;
	t113 = sin(qJ(2));
	t116 = t112 * t113;
	t114 = cos(qJ(2));
	t115 = t112 * t114;
	t1 = [t120, t109 * t116 - t111 * t114, t109 * t115 + t111 * t113, (t111 * pkin(2) + qJ(3) * t119) * t114 + (-pkin(2) * t119 + t111 * qJ(3)) * t113 + pkin(7) * t120 + t111 * pkin(1) + 0; -t118, -t109 * t114 - t111 * t116, t109 * t113 - t111 * t115, (t109 * pkin(2) - qJ(3) * t117) * t114 + (pkin(2) * t117 + t109 * qJ(3)) * t113 - pkin(7) * t118 + t109 * pkin(1) + 0; t112, -t110 * t113, -t110 * t114, t112 * pkin(7) + qJ(1) + 0 + (pkin(2) * t113 - qJ(3) * t114) * t110; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:02
	% EndTime: 2020-11-04 21:01:02
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (36->28), mult. (64->48), div. (0->0), fcn. (85->8), ass. (0->22)
	t123 = sin(pkin(10));
	t141 = t123 * qJ(3);
	t124 = sin(pkin(6));
	t127 = sin(qJ(4));
	t140 = t124 * t127;
	t129 = cos(qJ(4));
	t139 = t124 * t129;
	t130 = cos(qJ(2));
	t138 = t124 * t130;
	t131 = pkin(3) + pkin(7);
	t137 = t124 * t131;
	t125 = cos(pkin(10));
	t136 = t125 * qJ(3);
	t126 = cos(pkin(6));
	t128 = sin(qJ(2));
	t135 = t126 * t128;
	t134 = t126 * t130;
	t132 = pkin(2) + pkin(8);
	t133 = t126 * t132;
	t122 = t123 * t134 + t125 * t128;
	t121 = t123 * t128 - t125 * t134;
	t1 = [t122 * t127 + t123 * t139, t122 * t129 - t123 * t140, -t123 * t135 + t125 * t130, (t125 * t132 + t126 * t141) * t130 + (-t123 * t133 + t136) * t128 + t123 * t137 + t125 * pkin(1) + 0; t121 * t127 - t125 * t139, t121 * t129 + t125 * t140, t123 * t130 + t125 * t135, (t123 * t132 - t126 * t136) * t130 + (t125 * t133 + t141) * t128 - t125 * t137 + t123 * pkin(1) + 0; t126 * t129 - t127 * t138, -t126 * t127 - t129 * t138, t124 * t128, t131 * t126 + qJ(1) + 0 + (-qJ(3) * t130 + t128 * t132) * t124; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:02
	% EndTime: 2020-11-04 21:01:02
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (52->27), mult. (102->41), div. (0->0), fcn. (123->8), ass. (0->21)
	t145 = sin(pkin(6));
	t147 = cos(pkin(6));
	t149 = sin(qJ(2));
	t151 = cos(qJ(2));
	t153 = pkin(2) + pkin(8);
	t148 = sin(qJ(4));
	t150 = cos(qJ(4));
	t157 = t148 * pkin(4) - qJ(5) * t150 + qJ(3);
	t155 = -t153 * t149 + t157 * t151;
	t156 = t150 * pkin(4) + qJ(5) * t148 + pkin(3) + pkin(7);
	t165 = t156 * t145 + t155 * t147;
	t162 = t147 * t149;
	t161 = t147 * t151;
	t160 = t148 * t145;
	t159 = t150 * t145;
	t154 = t157 * t149 + t153 * t151 + pkin(1);
	t146 = cos(pkin(10));
	t144 = sin(pkin(10));
	t143 = t144 * t161 + t146 * t149;
	t142 = -t144 * t149 + t146 * t161;
	t1 = [-t144 * t162 + t146 * t151, -t143 * t148 - t144 * t159, -t143 * t150 + t144 * t160, t165 * t144 + t154 * t146 + 0; t144 * t151 + t146 * t162, t142 * t148 + t146 * t159, t142 * t150 - t146 * t160, t154 * t144 - t165 * t146 + 0; t145 * t149, -t147 * t150 + t151 * t160, t147 * t148 + t151 * t159, -t155 * t145 + t156 * t147 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:02
	% EndTime: 2020-11-04 21:01:02
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (78->35), mult. (140->56), div. (0->0), fcn. (174->10), ass. (0->31)
	t173 = sin(pkin(6));
	t175 = cos(pkin(6));
	t177 = sin(qJ(4));
	t180 = cos(qJ(4));
	t182 = pkin(4) + pkin(9);
	t186 = qJ(5) * t177 + t182 * t180 + pkin(3) + pkin(7);
	t171 = pkin(2) + pkin(5) + pkin(8);
	t178 = sin(qJ(2));
	t181 = cos(qJ(2));
	t187 = qJ(5) * t180 - t182 * t177 - qJ(3);
	t197 = t171 * t178 + t187 * t181;
	t199 = -t186 * t173 + t197 * t175;
	t195 = t173 * t178;
	t194 = t173 * t180;
	t193 = t175 * t178;
	t192 = t175 * t181;
	t191 = t177 * t173;
	t190 = t178 * t180;
	t189 = t180 * t181;
	t188 = t175 * t189 - t191;
	t184 = t171 * t181 - t187 * t178 + pkin(1);
	t179 = cos(qJ(6));
	t176 = sin(qJ(6));
	t174 = cos(pkin(10));
	t172 = sin(pkin(10));
	t170 = t173 * t189 + t175 * t177;
	t169 = t172 * t193 - t174 * t181;
	t168 = t172 * t181 + t174 * t193;
	t167 = -t172 * t190 + t188 * t174;
	t166 = -t188 * t172 - t174 * t190;
	t1 = [t166 * t176 - t179 * t169, t166 * t179 + t176 * t169, t172 * t194 + (t172 * t192 + t174 * t178) * t177, -t199 * t172 + t184 * t174 + 0; t167 * t176 + t168 * t179, t167 * t179 - t168 * t176, -t174 * t194 - (-t172 * t178 + t174 * t192) * t177, t184 * t172 + t199 * t174 + 0; t170 * t176 + t179 * t195, t170 * t179 - t176 * t195, t175 * t180 - t181 * t191, t173 * t197 + t186 * t175 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end