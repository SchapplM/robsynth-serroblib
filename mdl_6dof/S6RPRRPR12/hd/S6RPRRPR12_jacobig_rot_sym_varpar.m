% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:43
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRRPR12_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRRPR12_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:25
	% EndTime: 2019-10-10 01:43:25
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
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
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
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->13), mult. (46->31), div. (0->0), fcn. (67->10), ass. (0->19)
	t141 = sin(pkin(6));
	t146 = sin(qJ(1));
	t154 = t141 * t146;
	t148 = cos(qJ(1));
	t153 = t141 * t148;
	t139 = sin(pkin(12));
	t152 = t146 * t139;
	t142 = cos(pkin(12));
	t151 = t146 * t142;
	t150 = t148 * t139;
	t149 = t148 * t142;
	t147 = cos(qJ(3));
	t145 = sin(qJ(3));
	t144 = cos(pkin(6));
	t143 = cos(pkin(7));
	t140 = sin(pkin(7));
	t138 = -t144 * t151 - t150;
	t137 = t144 * t149 - t152;
	t1 = [0, 0, -t138 * t140 + t143 * t154, (-t144 * t152 + t149) * t145 + (-t138 * t143 - t140 * t154) * t147, 0, 0; 0, 0, -t137 * t140 - t143 * t153, (t144 * t150 + t151) * t145 + (-t137 * t143 + t140 * t153) * t147, 0, 0; 1, 0, -t141 * t142 * t140 + t144 * t143, -t144 * t140 * t147 + (-t142 * t143 * t147 + t139 * t145) * t141, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:43:26
	% EndTime: 2019-10-10 01:43:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->20), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->30)
	t171 = sin(pkin(7));
	t175 = cos(pkin(6));
	t191 = t171 * t175;
	t172 = sin(pkin(6));
	t178 = sin(qJ(1));
	t190 = t172 * t178;
	t181 = cos(qJ(1));
	t189 = t172 * t181;
	t173 = cos(pkin(12));
	t174 = cos(pkin(7));
	t188 = t173 * t174;
	t170 = sin(pkin(12));
	t187 = t178 * t170;
	t186 = t178 * t173;
	t185 = t181 * t170;
	t184 = t181 * t173;
	t166 = t175 * t184 - t187;
	t183 = -t166 * t174 + t171 * t189;
	t168 = -t175 * t186 - t185;
	t182 = t168 * t174 + t171 * t190;
	t180 = cos(qJ(3));
	t179 = cos(qJ(4));
	t177 = sin(qJ(3));
	t176 = sin(qJ(4));
	t169 = -t175 * t187 + t184;
	t167 = t175 * t185 + t186;
	t165 = -t172 * t173 * t171 + t175 * t174;
	t164 = -t168 * t171 + t174 * t190;
	t163 = -t166 * t171 - t174 * t189;
	t1 = [0, 0, t164, t169 * t177 - t182 * t180, 0, (t169 * t180 + t182 * t177) * t179 + t164 * t176; 0, 0, t163, t167 * t177 + t183 * t180, 0, (t167 * t180 - t183 * t177) * t179 + t163 * t176; 1, 0, t165, -t180 * t191 + (t170 * t177 - t180 * t188) * t172, 0, (t177 * t191 + (t170 * t180 + t177 * t188) * t172) * t179 + t165 * t176;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end