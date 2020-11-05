% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR3 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:02
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:54
	% EndTime: 2020-11-04 22:02:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:54
	% EndTime: 2020-11-04 22:02:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t99 = cos(qJ(1));
	t98 = sin(qJ(1));
	t1 = [t99, -t98, 0, 0; t98, t99, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:54
	% EndTime: 2020-11-04 22:02:54
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t100 = sin(pkin(6));
	t103 = sin(qJ(1));
	t111 = t103 * t100;
	t102 = sin(qJ(2));
	t110 = t103 * t102;
	t104 = cos(qJ(2));
	t109 = t103 * t104;
	t105 = cos(qJ(1));
	t108 = t105 * t100;
	t107 = t105 * t102;
	t106 = t105 * t104;
	t101 = cos(pkin(6));
	t1 = [-t101 * t110 + t106, -t101 * t109 - t107, t111, t105 * pkin(1) + pkin(8) * t111 + 0; t101 * t107 + t109, t101 * t106 - t110, -t108, t103 * pkin(1) - pkin(8) * t108 + 0; t100 * t102, t100 * t104, t101, t101 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:54
	% EndTime: 2020-11-04 22:02:54
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (33->28), div. (0->0), fcn. (44->12), ass. (0->20)
	t131 = pkin(2) * sin(qJ(2));
	t124 = qJ(2) + pkin(11);
	t130 = cos(qJ(1));
	t129 = sin(qJ(1));
	t127 = pkin(8) + qJ(3);
	t126 = cos(pkin(6));
	t125 = sin(pkin(6));
	t123 = cos(t124);
	t122 = sin(t124);
	t121 = pkin(6) - t124;
	t120 = pkin(6) + t124;
	t119 = cos(qJ(2)) * pkin(2) + pkin(1);
	t118 = cos(t121);
	t117 = cos(t120);
	t116 = sin(t121);
	t115 = sin(t120);
	t114 = t117 + t118;
	t113 = -t115 + t116;
	t112 = -t125 * t127 + t126 * t131;
	t1 = [t130 * t123 + t129 * t113 / 0.2e1, -t130 * t122 - t129 * t114 / 0.2e1, t129 * t125, -t112 * t129 + t130 * t119 + 0; t129 * t123 - t130 * t113 / 0.2e1, -t129 * t122 + t130 * t114 / 0.2e1, -t130 * t125, t130 * t112 + t129 * t119 + 0; t118 / 0.2e1 - t117 / 0.2e1, t116 / 0.2e1 + t115 / 0.2e1, t126, t125 * t131 + t126 * t127 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:54
	% EndTime: 2020-11-04 22:02:54
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->30), mult. (81->43), div. (0->0), fcn. (101->16), ass. (0->30)
	t143 = qJ(2) + pkin(11);
	t141 = sin(t143);
	t146 = sin(pkin(6));
	t160 = t141 * t146;
	t152 = sin(qJ(1));
	t159 = t146 * t152;
	t154 = cos(qJ(1));
	t158 = t146 * t154;
	t157 = t152 * t141;
	t156 = t154 * t141;
	t145 = sin(pkin(11));
	t148 = cos(pkin(11));
	t137 = t148 * pkin(3) + qJ(4) * t145 + pkin(2);
	t138 = -t145 * pkin(3) + qJ(4) * t148;
	t151 = sin(qJ(2));
	t153 = cos(qJ(2));
	t155 = t137 * t151 - t138 * t153;
	t150 = pkin(8) + qJ(3);
	t149 = cos(pkin(6));
	t147 = cos(pkin(12));
	t144 = sin(pkin(12));
	t142 = cos(t143);
	t140 = pkin(6) - t143;
	t139 = pkin(6) + t143;
	t136 = cos(t139) + cos(t140);
	t135 = t152 * t142 + t149 * t156;
	t134 = -t154 * t142 + t149 * t157;
	t133 = t137 * t153 + t138 * t151 + pkin(1);
	t132 = t146 * t150 - t155 * t149;
	t1 = [-t134 * t147 + t144 * t159, t134 * t144 + t147 * t159, t156 + t152 * t136 / 0.2e1, t132 * t152 + t133 * t154 + 0; t135 * t147 - t144 * t158, -t135 * t144 - t147 * t158, t157 - t154 * t136 / 0.2e1, -t132 * t154 + t133 * t152 + 0; t149 * t144 + t147 * t160, -t144 * t160 + t149 * t147, -sin(t140) / 0.2e1 - sin(t139) / 0.2e1, t155 * t146 + t149 * t150 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:54
	% EndTime: 2020-11-04 22:02:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (108->34), mult. (94->45), div. (0->0), fcn. (114->18), ass. (0->33)
	t177 = qJ(2) + pkin(11);
	t173 = sin(t177);
	t179 = sin(pkin(6));
	t192 = t173 * t179;
	t184 = sin(qJ(1));
	t191 = t179 * t184;
	t186 = cos(qJ(1));
	t190 = t179 * t186;
	t189 = t184 * t173;
	t188 = t186 * t173;
	t169 = cos(pkin(12)) * pkin(4) + pkin(3);
	t178 = sin(pkin(11));
	t180 = cos(pkin(11));
	t182 = pkin(9) + qJ(4);
	t165 = t169 * t180 + t182 * t178 + pkin(2);
	t166 = -t178 * t169 + t182 * t180;
	t183 = sin(qJ(2));
	t185 = cos(qJ(2));
	t187 = t165 * t183 - t166 * t185;
	t181 = cos(pkin(6));
	t176 = pkin(12) + qJ(5);
	t175 = cos(t177);
	t174 = cos(t176);
	t172 = sin(t176);
	t171 = pkin(6) - t177;
	t170 = pkin(6) + t177;
	t168 = sin(pkin(12)) * pkin(4) + qJ(3) + pkin(8);
	t167 = cos(t170) + cos(t171);
	t164 = t184 * t175 + t181 * t188;
	t163 = -t186 * t175 + t181 * t189;
	t162 = t165 * t185 + t166 * t183 + pkin(1);
	t161 = t179 * t168 - t187 * t181;
	t1 = [-t163 * t174 + t172 * t191, t163 * t172 + t174 * t191, t188 + t184 * t167 / 0.2e1, t161 * t184 + t162 * t186 + 0; t164 * t174 - t172 * t190, -t164 * t172 - t174 * t190, t189 - t186 * t167 / 0.2e1, -t161 * t186 + t162 * t184 + 0; t181 * t172 + t174 * t192, -t172 * t192 + t181 * t174, -sin(t171) / 0.2e1 - sin(t170) / 0.2e1, t168 * t181 + t187 * t179 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:02:54
	% EndTime: 2020-11-04 22:02:54
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (187->47), mult. (165->63), div. (0->0), fcn. (207->20), ass. (0->44)
	t218 = qJ(2) + pkin(11);
	t214 = sin(t218);
	t220 = sin(pkin(6));
	t235 = t214 * t220;
	t226 = sin(qJ(1));
	t234 = t220 * t226;
	t229 = cos(qJ(1));
	t233 = t220 * t229;
	t232 = t226 * t214;
	t231 = t229 * t214;
	t210 = cos(pkin(12)) * pkin(4) + pkin(3);
	t219 = sin(pkin(11));
	t221 = cos(pkin(11));
	t223 = pkin(9) + qJ(4);
	t205 = t210 * t221 + t223 * t219 + pkin(2);
	t206 = -t219 * t210 + t223 * t221;
	t225 = sin(qJ(2));
	t228 = cos(qJ(2));
	t230 = t205 * t225 - t206 * t228;
	t227 = cos(qJ(6));
	t224 = sin(qJ(6));
	t222 = cos(pkin(6));
	t217 = pkin(12) + qJ(5);
	t216 = cos(t218);
	t215 = cos(t217);
	t213 = sin(t217);
	t212 = pkin(6) - t218;
	t211 = pkin(6) + t218;
	t209 = sin(pkin(12)) * pkin(4) + qJ(3) + pkin(8);
	t208 = cos(t211) + cos(t212);
	t207 = -sin(t212) / 0.2e1 - sin(t211) / 0.2e1;
	t204 = t226 * t216 + t222 * t231;
	t203 = -t229 * t216 + t222 * t232;
	t202 = t222 * t213 + t215 * t235;
	t201 = t213 * t235 - t222 * t215;
	t200 = t231 + t226 * t208 / 0.2e1;
	t199 = t232 - t229 * t208 / 0.2e1;
	t198 = -t203 * t215 + t213 * t234;
	t197 = t204 * t215 - t213 * t233;
	t196 = t204 * t213 + t215 * t233;
	t195 = t203 * t213 + t215 * t234;
	t194 = t205 * t228 + t206 * t225 + pkin(1);
	t193 = t220 * t209 - t230 * t222;
	t1 = [t198 * t227 + t200 * t224, -t198 * t224 + t200 * t227, -t195, t198 * pkin(5) - t195 * pkin(10) + t193 * t226 + t194 * t229 + 0; t197 * t227 + t199 * t224, -t197 * t224 + t199 * t227, t196, t197 * pkin(5) + t196 * pkin(10) - t193 * t229 + t194 * t226 + 0; t202 * t227 + t207 * t224, -t202 * t224 + t207 * t227, t201, t202 * pkin(5) + t201 * pkin(10) + t209 * t222 + t230 * t220 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end