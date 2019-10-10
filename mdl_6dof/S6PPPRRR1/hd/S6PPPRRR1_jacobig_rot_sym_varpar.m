% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:49
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPPRRR1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPPRRR1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobig_rot_sym_varpar: pkin has to be [14x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:21
	% EndTime: 2019-10-10 08:49:21
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (18->16), mult. (52->38), div. (0->0), fcn. (73->12), ass. (0->19)
	t67 = sin(pkin(12));
	t76 = cos(pkin(6));
	t80 = t67 * t76;
	t69 = sin(pkin(7));
	t70 = sin(pkin(6));
	t79 = t69 * t70;
	t75 = cos(pkin(7));
	t78 = t70 * t75;
	t73 = cos(pkin(12));
	t77 = t73 * t76;
	t74 = cos(pkin(8));
	t72 = cos(pkin(13));
	t71 = cos(pkin(14));
	t68 = sin(pkin(8));
	t66 = sin(pkin(13));
	t65 = sin(pkin(14));
	t64 = -t73 * t66 - t72 * t80;
	t63 = -t67 * t66 + t72 * t77;
	t1 = [0, 0, 0, -(-(-t66 * t80 + t73 * t72) * t65 + (t64 * t75 + t67 * t79) * t71) * t68 + (-t64 * t69 + t67 * t78) * t74, 0, 0; 0, 0, 0, -(-(t66 * t77 + t67 * t72) * t65 + (t63 * t75 - t73 * t79) * t71) * t68 + (-t63 * t69 - t73 * t78) * t74, 0, 0; 0, 0, 0, -(t76 * t69 * t71 + (t71 * t72 * t75 - t65 * t66) * t70) * t68 + (-t72 * t79 + t76 * t75) * t74, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:22
	% EndTime: 2019-10-10 08:49:22
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (49->26), mult. (144->58), div. (0->0), fcn. (199->14), ass. (0->33)
	t133 = sin(pkin(12));
	t142 = cos(pkin(6));
	t152 = t133 * t142;
	t135 = sin(pkin(7));
	t136 = sin(pkin(6));
	t151 = t135 * t136;
	t150 = t135 * t142;
	t141 = cos(pkin(7));
	t149 = t136 * t141;
	t138 = cos(pkin(13));
	t148 = t138 * t141;
	t139 = cos(pkin(12));
	t147 = t139 * t142;
	t132 = sin(pkin(13));
	t127 = -t133 * t132 + t138 * t147;
	t146 = t127 * t141 - t139 * t151;
	t129 = -t139 * t132 - t138 * t152;
	t145 = t129 * t141 + t133 * t151;
	t144 = cos(qJ(4));
	t143 = sin(qJ(4));
	t140 = cos(pkin(8));
	t137 = cos(pkin(14));
	t134 = sin(pkin(8));
	t131 = sin(pkin(14));
	t130 = -t132 * t152 + t139 * t138;
	t128 = t132 * t147 + t133 * t138;
	t126 = -t138 * t151 + t142 * t141;
	t125 = -t129 * t135 + t133 * t149;
	t124 = -t127 * t135 - t139 * t149;
	t123 = t137 * t150 + (-t131 * t132 + t137 * t148) * t136;
	t122 = -t130 * t131 + t145 * t137;
	t121 = -t128 * t131 + t146 * t137;
	t1 = [0, 0, 0, -t122 * t134 + t125 * t140, (t130 * t137 + t145 * t131) * t143 + (-t122 * t140 - t125 * t134) * t144, 0; 0, 0, 0, -t121 * t134 + t124 * t140, (t128 * t137 + t146 * t131) * t143 + (-t121 * t140 - t124 * t134) * t144, 0; 0, 0, 0, -t123 * t134 + t126 * t140, (t136 * t132 * t137 + (t136 * t148 + t150) * t131) * t143 + (-t123 * t140 - t126 * t134) * t144, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:49:22
	% EndTime: 2019-10-10 08:49:22
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (101->32), mult. (294->70), div. (0->0), fcn. (404->16), ass. (0->44)
	t186 = sin(pkin(12));
	t195 = cos(pkin(6));
	t210 = t186 * t195;
	t188 = sin(pkin(7));
	t189 = sin(pkin(6));
	t209 = t188 * t189;
	t208 = t188 * t195;
	t194 = cos(pkin(7));
	t207 = t189 * t194;
	t191 = cos(pkin(13));
	t206 = t191 * t194;
	t192 = cos(pkin(12));
	t205 = t192 * t195;
	t185 = sin(pkin(13));
	t181 = t185 * t205 + t186 * t191;
	t184 = sin(pkin(14));
	t190 = cos(pkin(14));
	t180 = -t185 * t186 + t191 * t205;
	t201 = t180 * t194 - t192 * t209;
	t170 = -t181 * t184 + t201 * t190;
	t177 = -t180 * t188 - t192 * t207;
	t187 = sin(pkin(8));
	t193 = cos(pkin(8));
	t204 = t170 * t193 + t177 * t187;
	t183 = -t185 * t210 + t191 * t192;
	t182 = -t185 * t192 - t191 * t210;
	t200 = t182 * t194 + t186 * t209;
	t172 = -t183 * t184 + t200 * t190;
	t178 = -t182 * t188 + t186 * t207;
	t203 = t172 * t193 + t178 * t187;
	t175 = t190 * t208 + (-t184 * t185 + t190 * t206) * t189;
	t179 = -t191 * t209 + t194 * t195;
	t202 = t175 * t193 + t179 * t187;
	t199 = cos(qJ(4));
	t198 = cos(qJ(5));
	t197 = sin(qJ(4));
	t196 = sin(qJ(5));
	t176 = t185 * t189 * t190 + (t189 * t206 + t208) * t184;
	t174 = -t175 * t187 + t179 * t193;
	t173 = t183 * t190 + t200 * t184;
	t171 = t181 * t190 + t201 * t184;
	t169 = -t172 * t187 + t178 * t193;
	t168 = -t170 * t187 + t177 * t193;
	t1 = [0, 0, 0, t169, t173 * t197 - t203 * t199, (t173 * t199 + t203 * t197) * t196 - t169 * t198; 0, 0, 0, t168, t171 * t197 - t204 * t199, (t171 * t199 + t204 * t197) * t196 - t168 * t198; 0, 0, 0, t174, t176 * t197 - t202 * t199, (t176 * t199 + t202 * t197) * t196 - t174 * t198;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end