% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRP1 (for one body)
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

function Tc_mdh = S6PRPRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_fkine_fixb_body_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');
Tc_mdh=NaN(4,4);
%% Symbolic Calculation
if link_index == 0
	% From fkine_0_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:23
	% EndTime: 2020-11-04 21:01:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 1
	% From fkine_1_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:23
	% EndTime: 2020-11-04 21:01:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (2->2), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t100 = cos(pkin(10));
	t99 = sin(pkin(10));
	t1 = [t100, -t99, 0, 0; t99, t100, 0, 0; 0, 0, 1, qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 2
	% From fkine_2_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:23
	% EndTime: 2020-11-04 21:01:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (11->11), mult. (23->19), div. (0->0), fcn. (36->6), ass. (0->11)
	t101 = sin(pkin(10));
	t102 = sin(pkin(6));
	t110 = t101 * t102;
	t103 = cos(pkin(10));
	t109 = t103 * t102;
	t104 = cos(pkin(6));
	t105 = sin(qJ(2));
	t108 = t104 * t105;
	t106 = cos(qJ(2));
	t107 = t104 * t106;
	t1 = [-t101 * t108 + t103 * t106, -t101 * t107 - t103 * t105, t110, t103 * pkin(1) + pkin(7) * t110 + 0; t101 * t106 + t103 * t108, -t101 * t105 + t103 * t107, -t109, t101 * pkin(1) - pkin(7) * t109 + 0; t102 * t105, t102 * t106, t104, t104 * pkin(7) + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 3
	% From fkine_3_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:24
	% EndTime: 2020-11-04 21:01:24
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (53->22), mult. (37->26), div. (0->0), fcn. (44->12), ass. (0->20)
	t130 = pkin(2) * sin(qJ(2));
	t123 = qJ(2) + pkin(11);
	t128 = pkin(7) + qJ(3);
	t127 = cos(pkin(6));
	t126 = cos(pkin(10));
	t125 = sin(pkin(6));
	t124 = sin(pkin(10));
	t122 = cos(t123);
	t121 = sin(t123);
	t120 = pkin(6) - t123;
	t119 = pkin(6) + t123;
	t118 = cos(qJ(2)) * pkin(2) + pkin(1);
	t117 = cos(t119);
	t116 = sin(t120);
	t115 = cos(t120) / 0.2e1;
	t114 = sin(t119) / 0.2e1;
	t113 = -t125 * t128 + t127 * t130;
	t112 = t117 / 0.2e1 + t115;
	t111 = t114 - t116 / 0.2e1;
	t1 = [-t124 * t111 + t126 * t122, -t124 * t112 - t126 * t121, t124 * t125, -t124 * t113 + t126 * t118 + 0; t126 * t111 + t124 * t122, t126 * t112 - t124 * t121, -t126 * t125, t126 * t113 + t124 * t118 + 0; t115 - t117 / 0.2e1, t116 / 0.2e1 + t114, t127, t125 * t130 + t127 * t128 + qJ(1) + 0; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 4
	% From fkine_4_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:24
	% EndTime: 2020-11-04 21:01:24
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (73->35), mult. (85->53), div. (0->0), fcn. (105->16), ass. (0->27)
	t142 = sin(pkin(10));
	t146 = cos(pkin(6));
	t156 = t142 * t146;
	t143 = sin(pkin(6));
	t147 = pkin(7) + qJ(3);
	t155 = t143 * t147;
	t148 = sin(qJ(4));
	t154 = t143 * t148;
	t150 = cos(qJ(4));
	t153 = t143 * t150;
	t145 = cos(pkin(10));
	t152 = t145 * t146;
	t140 = qJ(2) + pkin(11);
	t151 = cos(qJ(2));
	t149 = sin(qJ(2));
	t144 = cos(pkin(11));
	t141 = sin(pkin(11));
	t139 = cos(t140);
	t138 = sin(t140);
	t137 = pkin(6) - t140;
	t136 = pkin(6) + t140;
	t135 = -t141 * pkin(3) + t144 * pkin(8);
	t134 = t144 * pkin(3) + t141 * pkin(8) + pkin(2);
	t133 = cos(t136) + cos(t137);
	t132 = t138 * t152 + t142 * t139;
	t131 = t138 * t156 - t145 * t139;
	t1 = [-t131 * t150 + t142 * t154, t131 * t148 + t142 * t153, t145 * t138 + t142 * t133 / 0.2e1, (t145 * t134 + t135 * t156) * t151 + (-t134 * t156 + t145 * t135) * t149 + t142 * t155 + t145 * pkin(1) + 0; t132 * t150 - t145 * t154, -t132 * t148 - t145 * t153, t142 * t138 - t145 * t133 / 0.2e1, (t142 * t134 - t135 * t152) * t151 + (t134 * t152 + t142 * t135) * t149 - t145 * t155 + t142 * pkin(1) + 0; t138 * t153 + t146 * t148, -t138 * t154 + t146 * t150, -sin(t137) / 0.2e1 - sin(t136) / 0.2e1, t146 * t147 + qJ(1) + 0 + (t134 * t149 - t135 * t151) * t143; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 5
	% From fkine_5_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:24
	% EndTime: 2020-11-04 21:01:24
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (134->48), mult. (156->71), div. (0->0), fcn. (198->18), ass. (0->38)
	t177 = sin(pkin(10));
	t181 = cos(pkin(6));
	t193 = t177 * t181;
	t178 = sin(pkin(6));
	t182 = pkin(7) + qJ(3);
	t192 = t178 * t182;
	t184 = sin(qJ(4));
	t191 = t178 * t184;
	t187 = cos(qJ(4));
	t190 = t178 * t187;
	t180 = cos(pkin(10));
	t189 = t180 * t181;
	t175 = qJ(2) + pkin(11);
	t188 = cos(qJ(2));
	t186 = cos(qJ(5));
	t185 = sin(qJ(2));
	t183 = sin(qJ(5));
	t179 = cos(pkin(11));
	t176 = sin(pkin(11));
	t174 = cos(t175);
	t173 = sin(t175);
	t172 = pkin(6) - t175;
	t171 = pkin(6) + t175;
	t170 = -t176 * pkin(3) + t179 * pkin(8);
	t169 = t179 * pkin(3) + t176 * pkin(8) + pkin(2);
	t168 = cos(t171) + cos(t172);
	t167 = -sin(t172) / 0.2e1 - sin(t171) / 0.2e1;
	t166 = t173 * t190 + t181 * t184;
	t165 = t173 * t191 - t181 * t187;
	t164 = t173 * t189 + t177 * t174;
	t163 = t173 * t193 - t180 * t174;
	t162 = t180 * t173 + t177 * t168 / 0.2e1;
	t161 = t177 * t173 - t180 * t168 / 0.2e1;
	t160 = -t163 * t187 + t177 * t191;
	t159 = t164 * t187 - t180 * t191;
	t158 = t164 * t184 + t180 * t190;
	t157 = t163 * t184 + t177 * t190;
	t1 = [t160 * t186 + t162 * t183, -t160 * t183 + t162 * t186, -t157, t160 * pkin(4) - t157 * pkin(9) + (t180 * t169 + t170 * t193) * t188 + (-t169 * t193 + t180 * t170) * t185 + t177 * t192 + t180 * pkin(1) + 0; t159 * t186 + t161 * t183, -t159 * t183 + t161 * t186, t158, t159 * pkin(4) + t158 * pkin(9) + (t177 * t169 - t170 * t189) * t188 + (t169 * t189 + t177 * t170) * t185 - t180 * t192 + t177 * pkin(1) + 0; t166 * t186 + t167 * t183, -t166 * t183 + t167 * t186, t165, t166 * pkin(4) + t165 * pkin(9) + t181 * t182 + qJ(1) + 0 + (t169 * t185 - t170 * t188) * t178; 0, 0, 0, 1;];
	Tc_mdh = t1;
elseif link_index == 6
	% From fkine_6_floatb_twist_rotmat_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-11-04 21:01:24
	% EndTime: 2020-11-04 21:01:24
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (162->53), mult. (173->75), div. (0->0), fcn. (216->18), ass. (0->43)
	t213 = qJ(2) + pkin(11);
	t209 = pkin(6) + t213;
	t210 = pkin(6) - t213;
	t205 = cos(t209) + cos(t210);
	t211 = sin(t213);
	t215 = sin(pkin(10));
	t218 = cos(pkin(10));
	t198 = t215 * t211 - t218 * t205 / 0.2e1;
	t222 = sin(qJ(5));
	t235 = t198 * t222;
	t199 = t218 * t211 + t215 * t205 / 0.2e1;
	t234 = t199 * t222;
	t204 = -sin(t210) / 0.2e1 - sin(t209) / 0.2e1;
	t233 = t204 * t222;
	t219 = cos(pkin(6));
	t232 = t215 * t219;
	t216 = sin(pkin(6));
	t221 = pkin(7) + qJ(3);
	t231 = t216 * t221;
	t223 = sin(qJ(4));
	t230 = t216 * t223;
	t226 = cos(qJ(4));
	t229 = t216 * t226;
	t228 = t218 * t219;
	t227 = cos(qJ(2));
	t225 = cos(qJ(5));
	t224 = sin(qJ(2));
	t220 = -qJ(6) - pkin(9);
	t217 = cos(pkin(11));
	t214 = sin(pkin(11));
	t212 = cos(t213);
	t208 = t225 * pkin(5) + pkin(4);
	t207 = -t214 * pkin(3) + t217 * pkin(8);
	t206 = t217 * pkin(3) + t214 * pkin(8) + pkin(2);
	t203 = t211 * t229 + t219 * t223;
	t202 = t211 * t230 - t219 * t226;
	t201 = t211 * t228 + t215 * t212;
	t200 = t211 * t232 - t218 * t212;
	t197 = -t200 * t226 + t215 * t230;
	t196 = t201 * t226 - t218 * t230;
	t195 = t201 * t223 + t218 * t229;
	t194 = t200 * t223 + t215 * t229;
	t1 = [t197 * t225 + t234, -t197 * t222 + t199 * t225, -t194, t197 * t208 + t194 * t220 + pkin(5) * t234 + (t218 * t206 + t207 * t232) * t227 + (-t206 * t232 + t218 * t207) * t224 + t215 * t231 + t218 * pkin(1) + 0; t196 * t225 + t235, -t196 * t222 + t198 * t225, t195, t196 * t208 - t195 * t220 + pkin(5) * t235 + (t215 * t206 - t207 * t228) * t227 + (t206 * t228 + t215 * t207) * t224 - t218 * t231 + t215 * pkin(1) + 0; t203 * t225 + t233, -t203 * t222 + t204 * t225, t202, pkin(5) * t233 - t202 * t220 + t203 * t208 + t219 * t221 + qJ(1) + 0 + (t206 * t224 - t207 * t227) * t216; 0, 0, 0, 1;];
	Tc_mdh = t1;
end