% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S6RRPRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:06
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobiaD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:08
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:10
	% DurationCPUTime: 0.89s
	% Computational Cost: add. (4795->72), mult. (2858->158), div. (686->14), fcn. (3330->7), ass. (0->75)
	t125 = sin(qJ(1));
	t119 = t125 ^ 2;
	t117 = qJ(2) + pkin(10) + qJ(4);
	t115 = sin(t117);
	t110 = t115 ^ 2;
	t116 = cos(t117);
	t113 = 0.1e1 / t116 ^ 2;
	t161 = t110 * t113;
	t105 = t119 * t161 + 0.1e1;
	t109 = t115 * t110;
	t111 = t116 ^ 2;
	t112 = 0.1e1 / t116;
	t118 = qJD(2) + qJD(4);
	t159 = t112 * t115;
	t133 = t118 * (t109 * t112 / t111 + t159);
	t126 = cos(qJ(1));
	t151 = qJD(1) * t126;
	t143 = t125 * t151;
	t166 = 0.1e1 / t105 ^ 2 * (t119 * t133 + t143 * t161);
	t175 = -0.2e1 * t166;
	t103 = 0.1e1 / t105;
	t141 = 0.1e1 + t161;
	t173 = t125 * t141;
	t98 = t103 * t173;
	t174 = t125 * t98 - 0.1e1;
	t155 = 0.1e1 / t125 * t126;
	t124 = t126 ^ 2;
	t172 = qJD(1) * (0.1e1 / t119 * t124 + 0.1e1) * t155;
	t153 = t125 * t115;
	t102 = atan2(-t153, -t116);
	t101 = cos(t102);
	t100 = sin(t102);
	t146 = t100 * t153;
	t97 = -t101 * t116 - t146;
	t94 = 0.1e1 / t97;
	t95 = 0.1e1 / t97 ^ 2;
	t171 = t103 - 0.1e1;
	t157 = t116 * t118;
	t137 = t115 * t124 * t157;
	t144 = t115 * t151;
	t156 = t118 * t125;
	t162 = t101 * t115;
	t145 = t113 * t156;
	t89 = (-(-t116 * t156 - t144) * t112 + t110 * t145) * t103;
	t84 = (-t125 * t89 + t118) * t162 + (-t144 + (t89 - t156) * t116) * t100;
	t169 = t84 * t94 * t95;
	t92 = t110 * t124 * t95 + 0.1e1;
	t170 = (t95 * t137 + (-t124 * t169 - t95 * t143) * t110) / t92 ^ 2;
	t90 = 0.1e1 / t92;
	t168 = t90 * t95;
	t165 = t118 * t98;
	t163 = t126 * t95;
	t160 = t110 * t125;
	t158 = t115 * t126;
	t121 = 0.1e1 / t125 ^ 2;
	t154 = t121 * t124;
	t152 = qJD(1) * t125;
	t150 = 0.2e1 * t169;
	t108 = t111 * t154 + 0.1e1;
	t149 = 0.2e1 / t108 ^ 2 * (-t111 * t172 - t121 * t137);
	t148 = t94 * t170;
	t147 = t90 * t157;
	t142 = 0.2e1 * t95 * t170;
	t140 = 0.1e1 + t154;
	t139 = t112 * t175;
	t138 = t101 * t103 * t110 * t112;
	t136 = t141 * t126;
	t135 = t140 * t115;
	t106 = 0.1e1 / t108;
	t88 = (t171 * t115 * t100 - t125 * t138) * t126;
	t87 = t115 * t149 * t155 + (qJD(1) * t135 - t155 * t157) * t106;
	t86 = -t174 * t162 + (-t125 + t98) * t116 * t100;
	t85 = t173 * t175 + (qJD(1) * t136 + 0.2e1 * t125 * t133) * t103;
	t82 = (-t94 * t90 * t152 + (-0.2e1 * t148 + (-t118 * t86 - t84) * t168) * t126) * t116 + (t86 * t126 * t142 + (-t126 * t118 * t94 - ((-t125 * t85 - t151 * t98) * t101 + (t174 * t89 + t156 - t165) * t100) * t95 * t158 + (t126 * t150 + t95 * t152) * t86 - ((t85 - t151) * t100 + (t89 * t98 + t118 + (-t89 - t165) * t125) * t101) * t116 * t163) * t90) * t115;
	t1 = [t139 * t158 + (t118 * t136 - t152 * t159) * t103, t85, 0, t85, 0, 0; (-t94 * t147 + (0.2e1 * t148 + (qJD(1) * t88 + t84) * t168) * t115) * t125 + (-t88 * t95 * t147 + (t88 * t142 + (t88 * t150 + ((0.2e1 * t115 * t166 + t157 + (-t112 * t89 * t160 - t157) * t103) * t100 + (t139 * t160 + t89 * t115 + (t109 * t145 - (t89 - 0.2e1 * t156) * t115) * t103) * t101) * t163) * t90) * t115 + (-t94 + (-(t119 - t124) * t138 + t171 * t146) * t95) * t115 * t90 * qJD(1)) * t126, t82, 0, t82, 0, 0; t140 * t116 * t149 + (0.2e1 * t116 * t172 + t118 * t135) * t106, t87, 0, t87, 0, 0;];
	JaD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobiaD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:06:09
	% EndTime: 2019-10-10 10:06:10
	% DurationCPUTime: 1.08s
	% Computational Cost: add. (5302->95), mult. (3810->207), div. (753->12), fcn. (4455->9), ass. (0->95)
	t176 = sin(qJ(1));
	t171 = qJ(2) + pkin(10) + qJ(4);
	t169 = sin(t171);
	t165 = 0.1e1 / t169 ^ 2;
	t170 = cos(t171);
	t168 = t170 ^ 2;
	t222 = t165 * t168;
	t197 = 0.1e1 + t222;
	t237 = t176 * t197;
	t173 = t176 ^ 2;
	t162 = t173 * t222 + 0.1e1;
	t160 = 0.1e1 / t162;
	t164 = 0.1e1 / t169;
	t178 = cos(qJ(1));
	t209 = qJD(1) * t178;
	t198 = t170 * t209;
	t172 = qJD(2) + qJD(4);
	t217 = t172 * t176;
	t200 = t165 * t217;
	t134 = ((t169 * t217 - t198) * t164 + t168 * t200) * t160;
	t236 = -t134 + t217;
	t193 = qJD(1) * t169 + qJD(6);
	t216 = t172 * t178;
	t235 = -t170 * t216 + t193 * t176;
	t213 = t176 * t170;
	t159 = atan2(-t213, t169);
	t150 = cos(t159);
	t149 = sin(t159);
	t202 = t149 * t213;
	t144 = t150 * t169 - t202;
	t141 = 0.1e1 / t144;
	t177 = cos(qJ(6));
	t212 = t176 * t177;
	t175 = sin(qJ(6));
	t214 = t175 * t178;
	t156 = t169 * t214 + t212;
	t152 = 0.1e1 / t156;
	t142 = 0.1e1 / t144 ^ 2;
	t153 = 0.1e1 / t156 ^ 2;
	t234 = t160 - 0.1e1;
	t226 = t150 * t170;
	t129 = (-t134 * t176 + t172) * t226 + (t236 * t169 - t198) * t149;
	t233 = t129 * t141 * t142;
	t194 = qJD(6) * t169 + qJD(1);
	t189 = t194 * t178;
	t139 = t175 * t189 + t235 * t177;
	t211 = t177 * t178;
	t215 = t175 * t176;
	t155 = -t169 * t211 + t215;
	t151 = t155 ^ 2;
	t148 = t151 * t153 + 0.1e1;
	t225 = t153 * t155;
	t140 = -t235 * t175 + t177 * t189;
	t230 = t140 * t152 * t153;
	t232 = (t139 * t225 - t151 * t230) / t148 ^ 2;
	t167 = t170 * t168;
	t223 = t164 * t170;
	t187 = t172 * (-t164 * t165 * t167 - t223);
	t220 = t168 * t176;
	t191 = t209 * t220;
	t231 = (t165 * t191 + t173 * t187) / t162 ^ 2;
	t229 = t142 * t170;
	t228 = t142 * t178;
	t227 = t149 * t176;
	t224 = t155 * t175;
	t174 = t178 ^ 2;
	t221 = t168 * t174;
	t219 = t169 * t172;
	t218 = t170 * t172;
	t210 = qJD(1) * t176;
	t137 = t142 * t221 + 0.1e1;
	t208 = 0.2e1 * (-t221 * t233 + (-t169 * t174 * t218 - t191) * t142) / t137 ^ 2;
	t207 = 0.2e1 * t233;
	t206 = 0.2e1 * t232;
	t205 = -0.2e1 * t231;
	t204 = t170 * t231;
	t203 = t170 * t228;
	t201 = t164 * t220;
	t196 = t170 * t208;
	t195 = 0.2e1 * t155 * t230;
	t192 = t150 * t160 * t164 * t168;
	t190 = t197 * t178;
	t188 = t152 * t177 + t153 * t224;
	t186 = t188 * t178;
	t158 = -t169 * t215 + t211;
	t157 = t169 * t212 + t214;
	t146 = 0.1e1 / t148;
	t145 = t160 * t237;
	t135 = 0.1e1 / t137;
	t133 = (t234 * t170 * t149 + t176 * t192) * t178;
	t131 = t169 * t227 + t226 + (-t149 * t169 - t150 * t213) * t145;
	t130 = t205 * t237 + (qJD(1) * t190 + 0.2e1 * t176 * t187) * t160;
	t127 = t170 * t186 * t206 + (t186 * t219 + (t188 * t210 + ((qJD(6) * t152 + t195) * t175 + (-t139 * t175 + (-qJD(6) * t155 + t140) * t177) * t153) * t178) * t170) * t146;
	t126 = (t131 * t229 + t141 * t169) * t178 * t208 + ((t141 * t210 + (t131 * t172 + t129) * t228) * t169 + (-t141 * t216 - (-t130 * t150 * t176 + t236 * t149 + (t134 * t227 - t149 * t172 - t150 * t209) * t145) * t203 + (t142 * t210 + t178 * t207) * t131 - ((-t130 + t209) * t149 + ((t145 * t176 - 0.1e1) * t172 + (-t145 + t176) * t134) * t150) * t169 * t228) * t170) * t135;
	t1 = [0.2e1 * t164 * t178 * t204 + (t172 * t190 + t210 * t223) * t160, t130, 0, t130, 0, 0; (t141 * t196 + (t141 * t219 + (qJD(1) * t133 + t129) * t229) * t135) * t176 + (t142 * t196 * t133 + (-((-0.2e1 * t204 + t219 + (-t134 * t201 - t219) * t160) * t149 + (t201 * t205 - t134 * t170 + (-t167 * t200 + (t134 - 0.2e1 * t217) * t170) * t160) * t150) * t203 + (t142 * t219 + t170 * t207) * t133 + (-t141 + ((t173 - t174) * t192 + t234 * t202) * t142) * t170 * qJD(1)) * t135) * t178, t126, 0, t126, 0, 0; (-t152 * t157 + t158 * t225) * t206 + (t158 * t195 + (-t158 * t139 - t157 * t140 + t194 * t155 * t212 - (-t172 * t213 - t193 * t178) * t224) * t153 + (t193 * t211 + (-t194 * t175 + t177 * t218) * t176) * t152) * t146, t127, 0, t127, 0, -0.2e1 * t232 + 0.2e1 * (t139 * t146 * t153 + (-t146 * t230 - t153 * t232) * t155) * t155;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,6);
end