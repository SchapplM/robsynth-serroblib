% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRPR11_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR11_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
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
	% StartTime: 2019-10-10 01:41:20
	% EndTime: 2019-10-10 01:41:20
	% DurationCPUTime: 0.07s
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
	% StartTime: 2019-10-10 01:41:21
	% EndTime: 2019-10-10 01:41:21
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->13), mult. (46->31), div. (0->0), fcn. (67->10), ass. (0->19)
	t151 = sin(pkin(6));
	t156 = sin(qJ(1));
	t164 = t151 * t156;
	t158 = cos(qJ(1));
	t163 = t151 * t158;
	t149 = sin(pkin(12));
	t162 = t156 * t149;
	t152 = cos(pkin(12));
	t161 = t156 * t152;
	t160 = t158 * t149;
	t159 = t158 * t152;
	t157 = cos(qJ(3));
	t155 = sin(qJ(3));
	t154 = cos(pkin(6));
	t153 = cos(pkin(7));
	t150 = sin(pkin(7));
	t148 = -t154 * t161 - t160;
	t147 = t154 * t159 - t162;
	t1 = [0, 0, -t148 * t150 + t153 * t164, (-t154 * t162 + t159) * t155 + (-t148 * t153 - t150 * t164) * t157, 0, 0; 0, 0, -t147 * t150 - t153 * t163, (t154 * t160 + t161) * t155 + (-t147 * t153 + t150 * t163) * t157, 0, 0; 1, 0, -t151 * t152 * t150 + t154 * t153, -t154 * t150 * t157 + (-t152 * t153 * t157 + t149 * t155) * t151, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:41:21
	% EndTime: 2019-10-10 01:41:21
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (33->20), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->30)
	t175 = sin(pkin(7));
	t179 = cos(pkin(6));
	t195 = t175 * t179;
	t176 = sin(pkin(6));
	t182 = sin(qJ(1));
	t194 = t176 * t182;
	t185 = cos(qJ(1));
	t193 = t176 * t185;
	t177 = cos(pkin(12));
	t178 = cos(pkin(7));
	t192 = t177 * t178;
	t174 = sin(pkin(12));
	t191 = t182 * t174;
	t190 = t182 * t177;
	t189 = t185 * t174;
	t188 = t185 * t177;
	t170 = t179 * t188 - t191;
	t187 = -t170 * t178 + t175 * t193;
	t172 = -t179 * t190 - t189;
	t186 = t172 * t178 + t175 * t194;
	t184 = cos(qJ(3));
	t183 = cos(qJ(4));
	t181 = sin(qJ(3));
	t180 = sin(qJ(4));
	t173 = -t179 * t191 + t188;
	t171 = t179 * t189 + t190;
	t169 = -t176 * t177 * t175 + t179 * t178;
	t168 = -t172 * t175 + t178 * t194;
	t167 = -t170 * t175 - t178 * t193;
	t1 = [0, 0, t168, t173 * t181 - t186 * t184, 0, (t173 * t184 + t186 * t181) * t180 - t168 * t183; 0, 0, t167, t171 * t181 + t187 * t184, 0, (t171 * t184 - t187 * t181) * t180 - t167 * t183; 1, 0, t169, -t184 * t195 + (t174 * t181 - t184 * t192) * t176, 0, (t181 * t195 + (t174 * t184 + t181 * t192) * t176) * t180 - t169 * t183;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end