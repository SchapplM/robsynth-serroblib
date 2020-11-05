% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRRR4 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 22:18
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6RRPRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [12x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:13
	% EndTime: 2020-11-04 22:18:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:13
	% EndTime: 2020-11-04 22:18:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t101 = cos(qJ(1));
	t100 = sin(qJ(1));
	t1 = [t101, -t100, 0, 0; t100, t101, 0, 0; 0, 0, 1, pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:13
	% EndTime: 2020-11-04 22:18:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->17), div. (0->0), fcn. (36->6), ass. (0->13)
	t102 = sin(pkin(6));
	t105 = sin(qJ(1));
	t113 = t105 * t102;
	t104 = sin(qJ(2));
	t112 = t105 * t104;
	t106 = cos(qJ(2));
	t111 = t105 * t106;
	t107 = cos(qJ(1));
	t110 = t107 * t102;
	t109 = t107 * t104;
	t108 = t107 * t106;
	t103 = cos(pkin(6));
	t1 = [-t103 * t112 + t108, -t103 * t111 - t109, t113, t107 * pkin(1) + pkin(8) * t113 + 0; t103 * t109 + t111, t103 * t108 - t112, -t110, t105 * pkin(1) - pkin(8) * t110 + 0; t102 * t104, t102 * t106, t103, t103 * pkin(8) + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:13
	% EndTime: 2020-11-04 22:18:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (33->28), div. (0->0), fcn. (44->12), ass. (0->20)
	t133 = pkin(2) * sin(qJ(2));
	t126 = qJ(2) + pkin(12);
	t132 = cos(qJ(1));
	t131 = sin(qJ(1));
	t129 = pkin(8) + qJ(3);
	t128 = cos(pkin(6));
	t127 = sin(pkin(6));
	t125 = cos(t126);
	t124 = sin(t126);
	t123 = pkin(6) - t126;
	t122 = pkin(6) + t126;
	t121 = cos(qJ(2)) * pkin(2) + pkin(1);
	t120 = cos(t123);
	t119 = cos(t122);
	t118 = sin(t123);
	t117 = sin(t122);
	t116 = t119 + t120;
	t115 = -t117 + t118;
	t114 = -t127 * t129 + t128 * t133;
	t1 = [t132 * t125 + t131 * t115 / 0.2e1, -t132 * t124 - t131 * t116 / 0.2e1, t131 * t127, -t114 * t131 + t132 * t121 + 0; t131 * t125 - t132 * t115 / 0.2e1, -t131 * t124 + t132 * t116 / 0.2e1, -t132 * t127, t132 * t114 + t131 * t121 + 0; t120 / 0.2e1 - t119 / 0.2e1, t118 / 0.2e1 + t117 / 0.2e1, t128, t127 * t133 + t128 * t129 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:13
	% EndTime: 2020-11-04 22:18:13
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (73->30), mult. (81->43), div. (0->0), fcn. (101->16), ass. (0->30)
	t147 = sin(pkin(6));
	t151 = sin(qJ(4));
	t162 = t147 * t151;
	t154 = cos(qJ(4));
	t161 = t147 * t154;
	t156 = cos(qJ(1));
	t160 = t147 * t156;
	t145 = qJ(2) + pkin(12);
	t143 = sin(t145);
	t153 = sin(qJ(1));
	t159 = t153 * t143;
	t158 = t156 * t143;
	t146 = sin(pkin(12));
	t148 = cos(pkin(12));
	t139 = t148 * pkin(3) + t146 * pkin(9) + pkin(2);
	t140 = -t146 * pkin(3) + t148 * pkin(9);
	t152 = sin(qJ(2));
	t155 = cos(qJ(2));
	t157 = t139 * t152 - t140 * t155;
	t150 = pkin(8) + qJ(3);
	t149 = cos(pkin(6));
	t144 = cos(t145);
	t142 = pkin(6) - t145;
	t141 = pkin(6) + t145;
	t138 = cos(t141) + cos(t142);
	t137 = t153 * t144 + t149 * t158;
	t136 = -t156 * t144 + t149 * t159;
	t135 = t139 * t155 + t140 * t152 + pkin(1);
	t134 = t147 * t150 - t157 * t149;
	t1 = [-t136 * t154 + t153 * t162, t136 * t151 + t153 * t161, t158 + t153 * t138 / 0.2e1, t134 * t153 + t135 * t156 + 0; t137 * t154 - t151 * t160, -t137 * t151 - t154 * t160, t159 - t156 * t138 / 0.2e1, -t134 * t156 + t135 * t153 + 0; t143 * t161 + t149 * t151, -t143 * t162 + t149 * t154, -sin(t142) / 0.2e1 - sin(t141) / 0.2e1, t157 * t147 + t149 * t150 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:13
	% EndTime: 2020-11-04 22:18:13
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (108->34), mult. (105->45), div. (0->0), fcn. (125->18), ass. (0->33)
	t176 = qJ(2) + pkin(12);
	t172 = sin(t176);
	t179 = sin(pkin(6));
	t196 = t172 * t179;
	t184 = sin(qJ(1));
	t195 = t179 * t184;
	t187 = cos(qJ(1));
	t194 = t179 * t187;
	t193 = t184 * t172;
	t192 = t187 * t172;
	t191 = pkin(4) * cos(qJ(4)) + pkin(3);
	t190 = sin(qJ(4)) * pkin(4) + qJ(3) + pkin(8);
	t178 = sin(pkin(12));
	t180 = cos(pkin(12));
	t188 = pkin(10) + pkin(9);
	t167 = t188 * t178 + t191 * t180 + pkin(2);
	t168 = -t191 * t178 + t188 * t180;
	t183 = sin(qJ(2));
	t186 = cos(qJ(2));
	t189 = t167 * t183 - t168 * t186;
	t181 = cos(pkin(6));
	t177 = qJ(4) + qJ(5);
	t175 = cos(t177);
	t174 = sin(t177);
	t173 = cos(t176);
	t171 = pkin(6) - t176;
	t170 = pkin(6) + t176;
	t169 = cos(t170) + cos(t171);
	t166 = t184 * t173 + t181 * t192;
	t165 = -t187 * t173 + t181 * t193;
	t164 = t167 * t186 + t168 * t183 + pkin(1);
	t163 = t179 * t190 - t189 * t181;
	t1 = [-t165 * t175 + t174 * t195, t165 * t174 + t175 * t195, t192 + t184 * t169 / 0.2e1, t163 * t184 + t164 * t187 + 0; t166 * t175 - t174 * t194, -t166 * t174 - t175 * t194, t193 - t187 * t169 / 0.2e1, -t163 * t187 + t164 * t184 + 0; t181 * t174 + t175 * t196, -t174 * t196 + t181 * t175, -sin(t171) / 0.2e1 - sin(t170) / 0.2e1, t189 * t179 + t190 * t181 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 22:18:13
	% EndTime: 2020-11-04 22:18:14
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (187->47), mult. (176->63), div. (0->0), fcn. (218->20), ass. (0->44)
	t219 = qJ(2) + pkin(12);
	t215 = sin(t219);
	t222 = sin(pkin(6));
	t241 = t215 * t222;
	t228 = sin(qJ(1));
	t240 = t222 * t228;
	t232 = cos(qJ(1));
	t239 = t222 * t232;
	t238 = t228 * t215;
	t237 = t232 * t215;
	t236 = pkin(4) * cos(qJ(4)) + pkin(3);
	t235 = sin(qJ(4)) * pkin(4) + qJ(3) + pkin(8);
	t221 = sin(pkin(12));
	t223 = cos(pkin(12));
	t233 = pkin(10) + pkin(9);
	t209 = t233 * t221 + t236 * t223 + pkin(2);
	t210 = -t236 * t221 + t233 * t223;
	t227 = sin(qJ(2));
	t231 = cos(qJ(2));
	t234 = t209 * t227 - t210 * t231;
	t229 = cos(qJ(6));
	t225 = sin(qJ(6));
	t224 = cos(pkin(6));
	t220 = qJ(4) + qJ(5);
	t218 = cos(t220);
	t217 = sin(t220);
	t216 = cos(t219);
	t214 = pkin(6) - t219;
	t213 = pkin(6) + t219;
	t212 = cos(t213) + cos(t214);
	t211 = -sin(t214) / 0.2e1 - sin(t213) / 0.2e1;
	t208 = t228 * t216 + t224 * t237;
	t207 = -t232 * t216 + t224 * t238;
	t206 = t224 * t217 + t218 * t241;
	t205 = t217 * t241 - t224 * t218;
	t204 = t237 + t228 * t212 / 0.2e1;
	t203 = t238 - t232 * t212 / 0.2e1;
	t202 = -t207 * t218 + t217 * t240;
	t201 = t208 * t218 - t217 * t239;
	t200 = t208 * t217 + t218 * t239;
	t199 = t207 * t217 + t218 * t240;
	t198 = t209 * t231 + t210 * t227 + pkin(1);
	t197 = t222 * t235 - t234 * t224;
	t1 = [t202 * t229 + t204 * t225, -t202 * t225 + t204 * t229, -t199, t202 * pkin(5) - t199 * pkin(11) + t197 * t228 + t198 * t232 + 0; t201 * t229 + t203 * t225, -t201 * t225 + t203 * t229, t200, t201 * pkin(5) + t200 * pkin(11) - t197 * t232 + t198 * t228 + 0; t206 * t229 + t211 * t225, -t206 * t225 + t211 * t229, t205, t206 * pkin(5) + t205 * pkin(11) + t234 * t222 + t235 * t224 + pkin(7) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
end