% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:05
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRPRRR7_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR7_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobig_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(13)) * t18, 0, 0, 0, 0; 0, -cos(pkin(13)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t59 = sin(pkin(6));
	t1 = [0, sin(pkin(13)) * t59, 0, 0, 0, 0; 0, -cos(pkin(13)) * t59, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:02
	% EndTime: 2019-10-09 22:05:02
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (19->17), mult. (54->39), div. (0->0), fcn. (78->12), ass. (0->20)
	t116 = sin(pkin(13));
	t119 = sin(pkin(6));
	t131 = t116 * t119;
	t118 = sin(pkin(7));
	t130 = t118 * t119;
	t121 = cos(pkin(13));
	t129 = t121 * t119;
	t124 = cos(pkin(6));
	t125 = sin(qJ(2));
	t128 = t124 * t125;
	t126 = cos(qJ(2));
	t127 = t124 * t126;
	t123 = cos(pkin(7));
	t122 = cos(pkin(8));
	t120 = cos(pkin(14));
	t117 = sin(pkin(8));
	t115 = sin(pkin(14));
	t114 = -t116 * t127 - t121 * t125;
	t113 = -t116 * t125 + t121 * t127;
	t1 = [0, t131, 0, -(-(-t116 * t128 + t121 * t126) * t115 + (t114 * t123 + t116 * t130) * t120) * t117 + (-t114 * t118 + t123 * t131) * t122, 0, 0; 0, -t129, 0, -(-(t116 * t126 + t121 * t128) * t115 + (t113 * t123 - t118 * t129) * t120) * t117 + (-t113 * t118 - t123 * t129) * t122, 0, 0; 0, t124, 0, -(t124 * t118 * t120 + (t120 * t123 * t126 - t115 * t125) * t119) * t117 + (t124 * t123 - t126 * t130) * t122, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:03
	% EndTime: 2019-10-09 22:05:03
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (50->27), mult. (146->59), div. (0->0), fcn. (204->14), ass. (0->34)
	t180 = sin(pkin(13));
	t183 = sin(pkin(6));
	t201 = t180 * t183;
	t182 = sin(pkin(7));
	t200 = t182 * t183;
	t188 = cos(pkin(6));
	t199 = t182 * t188;
	t185 = cos(pkin(13));
	t198 = t185 * t183;
	t187 = cos(pkin(7));
	t192 = cos(qJ(2));
	t197 = t187 * t192;
	t190 = sin(qJ(2));
	t196 = t188 * t190;
	t195 = t188 * t192;
	t175 = -t180 * t190 + t185 * t195;
	t194 = t175 * t187 - t182 * t198;
	t177 = -t180 * t195 - t185 * t190;
	t193 = t177 * t187 + t180 * t200;
	t191 = cos(qJ(4));
	t189 = sin(qJ(4));
	t186 = cos(pkin(8));
	t184 = cos(pkin(14));
	t181 = sin(pkin(8));
	t179 = sin(pkin(14));
	t178 = -t180 * t196 + t185 * t192;
	t176 = t180 * t192 + t185 * t196;
	t174 = t188 * t187 - t192 * t200;
	t173 = -t177 * t182 + t187 * t201;
	t172 = -t175 * t182 - t187 * t198;
	t171 = t184 * t199 + (-t179 * t190 + t184 * t197) * t183;
	t170 = -t178 * t179 + t193 * t184;
	t169 = -t176 * t179 + t194 * t184;
	t1 = [0, t201, 0, -t170 * t181 + t173 * t186, (t178 * t184 + t193 * t179) * t189 + (-t170 * t186 - t173 * t181) * t191, 0; 0, -t198, 0, -t169 * t181 + t172 * t186, (t176 * t184 + t194 * t179) * t189 + (-t169 * t186 - t172 * t181) * t191, 0; 0, t188, 0, -t171 * t181 + t174 * t186, (t183 * t190 * t184 + (t183 * t197 + t199) * t179) * t189 + (-t171 * t186 - t174 * t181) * t191, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:05:04
	% EndTime: 2019-10-09 22:05:04
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (102->33), mult. (296->71), div. (0->0), fcn. (409->16), ass. (0->45)
	t231 = sin(pkin(13));
	t234 = sin(pkin(6));
	t257 = t231 * t234;
	t233 = sin(pkin(7));
	t256 = t233 * t234;
	t239 = cos(pkin(6));
	t255 = t233 * t239;
	t236 = cos(pkin(13));
	t254 = t236 * t234;
	t238 = cos(pkin(7));
	t245 = cos(qJ(2));
	t253 = t238 * t245;
	t242 = sin(qJ(2));
	t252 = t239 * t242;
	t251 = t239 * t245;
	t227 = t231 * t245 + t236 * t252;
	t230 = sin(pkin(14));
	t235 = cos(pkin(14));
	t226 = -t231 * t242 + t236 * t251;
	t247 = t226 * t238 - t233 * t254;
	t217 = -t227 * t230 + t247 * t235;
	t223 = -t226 * t233 - t238 * t254;
	t232 = sin(pkin(8));
	t237 = cos(pkin(8));
	t250 = t217 * t237 + t223 * t232;
	t229 = -t231 * t252 + t236 * t245;
	t228 = -t231 * t251 - t236 * t242;
	t246 = t228 * t238 + t231 * t256;
	t219 = -t229 * t230 + t246 * t235;
	t224 = -t228 * t233 + t238 * t257;
	t249 = t219 * t237 + t224 * t232;
	t221 = t235 * t255 + (-t230 * t242 + t235 * t253) * t234;
	t225 = t239 * t238 - t245 * t256;
	t248 = t221 * t237 + t225 * t232;
	t244 = cos(qJ(4));
	t243 = cos(qJ(5));
	t241 = sin(qJ(4));
	t240 = sin(qJ(5));
	t222 = t234 * t242 * t235 + (t234 * t253 + t255) * t230;
	t220 = t229 * t235 + t246 * t230;
	t218 = t227 * t235 + t247 * t230;
	t216 = -t221 * t232 + t225 * t237;
	t215 = -t219 * t232 + t224 * t237;
	t214 = -t217 * t232 + t223 * t237;
	t1 = [0, t257, 0, t215, t220 * t241 - t249 * t244, (t220 * t244 + t249 * t241) * t240 - t215 * t243; 0, -t254, 0, t214, t218 * t241 - t250 * t244, (t218 * t244 + t250 * t241) * t240 - t214 * t243; 0, t239, 0, t216, t222 * t241 - t248 * t244, (t222 * t244 + t248 * t241) * t240 - t216 * t243;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end