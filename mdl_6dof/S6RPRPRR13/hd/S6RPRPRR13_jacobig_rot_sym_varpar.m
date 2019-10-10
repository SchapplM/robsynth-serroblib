% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:07
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRPRR13_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR13_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
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
	% StartTime: 2019-10-10 01:07:44
	% EndTime: 2019-10-10 01:07:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t104 = sin(pkin(6));
	t106 = cos(pkin(7));
	t111 = t104 * t106;
	t105 = cos(pkin(12));
	t107 = cos(pkin(6));
	t110 = t105 * t107;
	t109 = cos(qJ(1));
	t108 = sin(qJ(1));
	t103 = sin(pkin(7));
	t102 = sin(pkin(12));
	t1 = [0, 0, -(-t109 * t102 - t108 * t110) * t103 + t108 * t111, 0, 0, 0; 0, 0, -(-t108 * t102 + t109 * t110) * t103 - t109 * t111, 0, 0, 0; 1, 0, -t104 * t105 * t103 + t107 * t106, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:45
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
	t1 = [0, 0, -t118 * t120 + t123 * t134, 0, (-t124 * t132 + t129) * t127 + (t118 * t123 + t120 * t134) * t125, 0; 0, 0, -t117 * t120 - t123 * t133, 0, (t124 * t130 + t131) * t127 + (t117 * t123 - t120 * t133) * t125, 0; 1, 0, -t121 * t122 * t120 + t124 * t123, 0, t124 * t120 * t125 + (t122 * t123 * t125 + t119 * t127) * t121, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:07:45
	% EndTime: 2019-10-10 01:07:45
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->20), mult. (98->45), div. (0->0), fcn. (140->12), ass. (0->30)
	t161 = sin(pkin(7));
	t165 = cos(pkin(6));
	t181 = t161 * t165;
	t162 = sin(pkin(6));
	t168 = sin(qJ(1));
	t180 = t162 * t168;
	t171 = cos(qJ(1));
	t179 = t162 * t171;
	t163 = cos(pkin(12));
	t164 = cos(pkin(7));
	t178 = t163 * t164;
	t160 = sin(pkin(12));
	t177 = t168 * t160;
	t176 = t168 * t163;
	t175 = t171 * t160;
	t174 = t171 * t163;
	t156 = t165 * t174 - t177;
	t173 = -t156 * t164 + t161 * t179;
	t158 = -t165 * t176 - t175;
	t172 = t158 * t164 + t161 * t180;
	t170 = cos(qJ(3));
	t169 = cos(qJ(5));
	t167 = sin(qJ(3));
	t166 = sin(qJ(5));
	t159 = -t165 * t177 + t174;
	t157 = t165 * t175 + t176;
	t155 = -t162 * t163 * t161 + t165 * t164;
	t154 = -t158 * t161 + t164 * t180;
	t153 = -t156 * t161 - t164 * t179;
	t1 = [0, 0, t154, 0, t159 * t170 + t172 * t167, t154 * t166 - (t159 * t167 - t172 * t170) * t169; 0, 0, t153, 0, t157 * t170 - t173 * t167, t153 * t166 - (t157 * t167 + t173 * t170) * t169; 1, 0, t155, 0, t167 * t181 + (t160 * t170 + t167 * t178) * t162, t155 * t166 - (-t170 * t181 + (t160 * t167 - t170 * t178) * t162) * t169;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end