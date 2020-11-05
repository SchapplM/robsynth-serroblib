% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRP2 (for one body)
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   index of the body frame to be returned (0=base).
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% Tc_mdh [4x4]
%   homogenous transformation matrices for body frame of "link_index"

% Quelle: HybrDyn-Toolbox
% Datum: 2020-11-04 21:01
% Revision: de51baf798caa2364afaf24686304d90a3288510 (2020-11-04)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Tc_mdh = S6PRPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:45
	% EndTime: 2020-11-04 21:01:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:45
	% EndTime: 2020-11-04 21:01:45
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t104 = cos(pkin(10));
	t103 = sin(pkin(10));
	t1 = [t104, -t103, 0, 0; t103, t104, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:45
	% EndTime: 2020-11-04 21:01:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t105 = sin(pkin(10));
	t106 = sin(pkin(6));
	t114 = t105 * t106;
	t107 = cos(pkin(10));
	t113 = t107 * t106;
	t108 = cos(pkin(6));
	t109 = sin(qJ(2));
	t112 = t108 * t109;
	t110 = cos(qJ(2));
	t111 = t108 * t110;
	t1 = [-t105 * t112 + t107 * t110, -t105 * t111 - t107 * t109, t114, t107 * pkin(1) + pkin(7) * t114 + 0; t105 * t110 + t107 * t112, -t105 * t109 + t107 * t111, -t113, t105 * pkin(1) - pkin(7) * t113 + 0; t106 * t109, t106 * t110, t108, t108 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:45
	% EndTime: 2020-11-04 21:01:45
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (37->26), div. (0->0), fcn. (44->12), ass. (0->20)
	t134 = pkin(2) * sin(qJ(2));
	t127 = qJ(2) + pkin(11);
	t132 = pkin(7) + qJ(3);
	t131 = cos(pkin(6));
	t130 = cos(pkin(10));
	t129 = sin(pkin(6));
	t128 = sin(pkin(10));
	t126 = cos(t127);
	t125 = sin(t127);
	t124 = pkin(6) - t127;
	t123 = pkin(6) + t127;
	t122 = cos(qJ(2)) * pkin(2) + pkin(1);
	t121 = cos(t123);
	t120 = sin(t124);
	t119 = cos(t124) / 0.2e1;
	t118 = sin(t123) / 0.2e1;
	t117 = -t129 * t132 + t131 * t134;
	t116 = t121 / 0.2e1 + t119;
	t115 = t118 - t120 / 0.2e1;
	t1 = [-t128 * t115 + t130 * t126, -t128 * t116 - t130 * t125, t128 * t129, -t128 * t117 + t130 * t122 + 0; t130 * t115 + t128 * t126, t130 * t116 - t128 * t125, -t130 * t129, t130 * t117 + t128 * t122 + 0; t119 - t121 / 0.2e1, t120 / 0.2e1 + t118, t131, t129 * t134 + t131 * t132 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:45
	% EndTime: 2020-11-04 21:01:45
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (73->35), mult. (85->53), div. (0->0), fcn. (105->16), ass. (0->27)
	t146 = sin(pkin(10));
	t150 = cos(pkin(6));
	t160 = t146 * t150;
	t147 = sin(pkin(6));
	t151 = pkin(7) + qJ(3);
	t159 = t147 * t151;
	t152 = sin(qJ(4));
	t158 = t147 * t152;
	t154 = cos(qJ(4));
	t157 = t147 * t154;
	t149 = cos(pkin(10));
	t156 = t149 * t150;
	t144 = qJ(2) + pkin(11);
	t155 = cos(qJ(2));
	t153 = sin(qJ(2));
	t148 = cos(pkin(11));
	t145 = sin(pkin(11));
	t143 = cos(t144);
	t142 = sin(t144);
	t141 = pkin(6) - t144;
	t140 = pkin(6) + t144;
	t139 = -t145 * pkin(3) + t148 * pkin(8);
	t138 = t148 * pkin(3) + t145 * pkin(8) + pkin(2);
	t137 = cos(t140) + cos(t141);
	t136 = t142 * t156 + t146 * t143;
	t135 = t142 * t160 - t149 * t143;
	t1 = [-t135 * t154 + t146 * t158, t135 * t152 + t146 * t157, t149 * t142 + t146 * t137 / 0.2e1, (t149 * t138 + t139 * t160) * t155 + (-t138 * t160 + t149 * t139) * t153 + t146 * t159 + t149 * pkin(1) + 0; t136 * t154 - t149 * t158, -t136 * t152 - t149 * t157, t146 * t142 - t149 * t137 / 0.2e1, (t146 * t138 - t139 * t156) * t155 + (t138 * t156 + t146 * t139) * t153 - t149 * t159 + t146 * pkin(1) + 0; t142 * t157 + t150 * t152, -t142 * t158 + t150 * t154, -sin(t141) / 0.2e1 - sin(t140) / 0.2e1, t150 * t151 + qJ(1) + 0 + (t138 * t153 - t139 * t155) * t147; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:45
	% EndTime: 2020-11-04 21:01:45
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (134->48), mult. (156->71), div. (0->0), fcn. (198->18), ass. (0->38)
	t181 = sin(pkin(10));
	t185 = cos(pkin(6));
	t197 = t181 * t185;
	t182 = sin(pkin(6));
	t186 = pkin(7) + qJ(3);
	t196 = t182 * t186;
	t188 = sin(qJ(4));
	t195 = t182 * t188;
	t191 = cos(qJ(4));
	t194 = t182 * t191;
	t184 = cos(pkin(10));
	t193 = t184 * t185;
	t179 = qJ(2) + pkin(11);
	t192 = cos(qJ(2));
	t190 = cos(qJ(5));
	t189 = sin(qJ(2));
	t187 = sin(qJ(5));
	t183 = cos(pkin(11));
	t180 = sin(pkin(11));
	t178 = cos(t179);
	t177 = sin(t179);
	t176 = pkin(6) - t179;
	t175 = pkin(6) + t179;
	t174 = -t180 * pkin(3) + t183 * pkin(8);
	t173 = t183 * pkin(3) + t180 * pkin(8) + pkin(2);
	t172 = cos(t175) + cos(t176);
	t171 = -sin(t176) / 0.2e1 - sin(t175) / 0.2e1;
	t170 = t177 * t194 + t185 * t188;
	t169 = t177 * t195 - t185 * t191;
	t168 = t177 * t193 + t181 * t178;
	t167 = t177 * t197 - t184 * t178;
	t166 = t184 * t177 + t181 * t172 / 0.2e1;
	t165 = t181 * t177 - t184 * t172 / 0.2e1;
	t164 = -t167 * t191 + t181 * t195;
	t163 = t168 * t191 - t184 * t195;
	t162 = t168 * t188 + t184 * t194;
	t161 = t167 * t188 + t181 * t194;
	t1 = [t164 * t190 + t166 * t187, -t164 * t187 + t166 * t190, -t161, t164 * pkin(4) - t161 * pkin(9) + (t184 * t173 + t174 * t197) * t192 + (-t173 * t197 + t184 * t174) * t189 + t181 * t196 + t184 * pkin(1) + 0; t163 * t190 + t165 * t187, -t163 * t187 + t165 * t190, t162, t163 * pkin(4) + t162 * pkin(9) + (t181 * t173 - t174 * t193) * t192 + (t173 * t193 + t181 * t174) * t189 - t184 * t196 + t181 * pkin(1) + 0; t170 * t190 + t171 * t187, -t170 * t187 + t171 * t190, t169, t170 * pkin(4) + t169 * pkin(9) + t185 * t186 + qJ(1) + 0 + (t173 * t189 - t174 * t192) * t182; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:45
	% EndTime: 2020-11-04 21:01:46
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (204->54), mult. (220->77), div. (0->0), fcn. (280->18), ass. (0->44)
	t224 = sin(pkin(10));
	t228 = cos(pkin(6));
	t240 = t224 * t228;
	t225 = sin(pkin(6));
	t229 = pkin(7) + qJ(3);
	t239 = t225 * t229;
	t231 = sin(qJ(4));
	t238 = t225 * t231;
	t234 = cos(qJ(4));
	t237 = t225 * t234;
	t227 = cos(pkin(10));
	t236 = t227 * t228;
	t222 = qJ(2) + pkin(11);
	t235 = cos(qJ(2));
	t233 = cos(qJ(5));
	t232 = sin(qJ(2));
	t230 = sin(qJ(5));
	t226 = cos(pkin(11));
	t223 = sin(pkin(11));
	t221 = cos(t222);
	t220 = sin(t222);
	t219 = pkin(6) - t222;
	t218 = pkin(6) + t222;
	t217 = -t223 * pkin(3) + t226 * pkin(8);
	t216 = t226 * pkin(3) + t223 * pkin(8) + pkin(2);
	t215 = cos(t218) + cos(t219);
	t214 = -sin(t219) / 0.2e1 - sin(t218) / 0.2e1;
	t213 = t220 * t237 + t228 * t231;
	t212 = t220 * t238 - t228 * t234;
	t211 = t220 * t236 + t224 * t221;
	t210 = t220 * t240 - t227 * t221;
	t209 = t227 * t220 + t224 * t215 / 0.2e1;
	t208 = t224 * t220 - t227 * t215 / 0.2e1;
	t207 = -t210 * t234 + t224 * t238;
	t206 = t211 * t234 - t227 * t238;
	t205 = t211 * t231 + t227 * t237;
	t204 = t210 * t231 + t224 * t237;
	t203 = t213 * t233 + t214 * t230;
	t202 = t213 * t230 - t214 * t233;
	t201 = t207 * t233 + t209 * t230;
	t200 = t207 * t230 - t209 * t233;
	t199 = t206 * t233 + t208 * t230;
	t198 = t206 * t230 - t208 * t233;
	t1 = [t201, -t204, t200, t201 * pkin(5) + t200 * qJ(6) + t207 * pkin(4) - t204 * pkin(9) + (t227 * t216 + t217 * t240) * t235 + (-t216 * t240 + t227 * t217) * t232 + t224 * t239 + t227 * pkin(1) + 0; t199, t205, t198, t199 * pkin(5) + t198 * qJ(6) + t206 * pkin(4) + t205 * pkin(9) + (t224 * t216 - t217 * t236) * t235 + (t216 * t236 + t224 * t217) * t232 - t227 * t239 + t224 * pkin(1) + 0; t203, t212, t202, t213 * pkin(4) + t203 * pkin(5) + t212 * pkin(9) + t202 * qJ(6) + t228 * t229 + qJ(1) + 0 + (t216 * t232 - t217 * t235) * t225; 0, 0, 0, 1;];
	Tc_mdh = t1;
end