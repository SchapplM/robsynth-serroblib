% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRRP12
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 08:58
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRRP12_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRRP12_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t78 = sin(pkin(6));
	t80 = cos(pkin(7));
	t85 = t78 * t80;
	t79 = cos(pkin(12));
	t81 = cos(pkin(6));
	t84 = t79 * t81;
	t83 = cos(qJ(1));
	t82 = sin(qJ(1));
	t77 = sin(pkin(7));
	t76 = sin(pkin(12));
	t1 = [0, 0, -(-t83 * t76 - t82 * t84) * t77 + t82 * t85, 0, 0, 0; 0, 0, -(-t82 * t76 + t83 * t84) * t77 - t83 * t85, 0, 0, 0; 1, 0, -t78 * t79 * t77 + t81 * t80, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:24
	% EndTime: 2019-10-10 08:58:24
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->13), mult. (46->31), div. (0->0), fcn. (67->10), ass. (0->19)
	t121 = sin(pkin(6));
	t126 = sin(qJ(1));
	t134 = t121 * t126;
	t128 = cos(qJ(1));
	t133 = t121 * t128;
	t119 = sin(pkin(12));
	t132 = t126 * t119;
	t122 = cos(pkin(12));
	t131 = t126 * t122;
	t130 = t128 * t119;
	t129 = t128 * t122;
	t127 = cos(qJ(3));
	t125 = sin(qJ(3));
	t124 = cos(pkin(6));
	t123 = cos(pkin(7));
	t120 = sin(pkin(7));
	t118 = -t124 * t131 - t130;
	t117 = t124 * t129 - t132;
	t1 = [0, 0, -t118 * t120 + t123 * t134, (-t124 * t132 + t129) * t125 + (-t118 * t123 - t120 * t134) * t127, 0, 0; 0, 0, -t117 * t120 - t123 * t133, (t124 * t130 + t131) * t125 + (-t117 * t123 + t120 * t133) * t127, 0, 0; 1, 0, -t121 * t122 * t120 + t124 * t123, -t124 * t120 * t127 + (-t122 * t123 * t127 + t119 * t125) * t121, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:25
	% EndTime: 2019-10-10 08:58:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->20), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->30)
	t163 = sin(pkin(7));
	t167 = cos(pkin(6));
	t183 = t163 * t167;
	t164 = sin(pkin(6));
	t170 = sin(qJ(1));
	t182 = t164 * t170;
	t173 = cos(qJ(1));
	t181 = t164 * t173;
	t165 = cos(pkin(12));
	t166 = cos(pkin(7));
	t180 = t165 * t166;
	t162 = sin(pkin(12));
	t179 = t170 * t162;
	t178 = t170 * t165;
	t177 = t173 * t162;
	t176 = t173 * t165;
	t158 = t167 * t176 - t179;
	t175 = -t158 * t166 + t163 * t181;
	t160 = -t167 * t178 - t177;
	t174 = t160 * t166 + t163 * t182;
	t172 = cos(qJ(3));
	t171 = cos(qJ(4));
	t169 = sin(qJ(3));
	t168 = sin(qJ(4));
	t161 = -t167 * t179 + t176;
	t159 = t167 * t177 + t178;
	t157 = -t164 * t165 * t163 + t167 * t166;
	t156 = -t160 * t163 + t166 * t182;
	t155 = -t158 * t163 - t166 * t181;
	t1 = [0, 0, t156, t161 * t169 - t174 * t172, (t161 * t172 + t174 * t169) * t168 - t156 * t171, 0; 0, 0, t155, t159 * t169 + t175 * t172, (t159 * t172 - t175 * t169) * t168 - t155 * t171, 0; 1, 0, t157, -t172 * t183 + (t162 * t169 - t172 * t180) * t164, (t169 * t183 + (t162 * t172 + t169 * t180) * t164) * t168 - t157 * t171, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 08:58:25
	% EndTime: 2019-10-10 08:58:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->20), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->30)
	t193 = sin(pkin(7));
	t197 = cos(pkin(6));
	t213 = t193 * t197;
	t194 = sin(pkin(6));
	t200 = sin(qJ(1));
	t212 = t194 * t200;
	t203 = cos(qJ(1));
	t211 = t194 * t203;
	t195 = cos(pkin(12));
	t196 = cos(pkin(7));
	t210 = t195 * t196;
	t192 = sin(pkin(12));
	t209 = t200 * t192;
	t208 = t200 * t195;
	t207 = t203 * t192;
	t206 = t203 * t195;
	t188 = t197 * t206 - t209;
	t205 = -t188 * t196 + t193 * t211;
	t190 = -t197 * t208 - t207;
	t204 = t190 * t196 + t193 * t212;
	t202 = cos(qJ(3));
	t201 = cos(qJ(4));
	t199 = sin(qJ(3));
	t198 = sin(qJ(4));
	t191 = -t197 * t209 + t206;
	t189 = t197 * t207 + t208;
	t187 = -t194 * t195 * t193 + t197 * t196;
	t186 = -t190 * t193 + t196 * t212;
	t185 = -t188 * t193 - t196 * t211;
	t1 = [0, 0, t186, t191 * t199 - t204 * t202, (t191 * t202 + t204 * t199) * t198 - t186 * t201, 0; 0, 0, t185, t189 * t199 + t205 * t202, (t189 * t202 - t205 * t199) * t198 - t185 * t201, 0; 1, 0, t187, -t202 * t213 + (t192 * t199 - t202 * t210) * t194, (t199 * t213 + (t192 * t202 + t199 * t210) * t194) * t198 - t187 * t201, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end