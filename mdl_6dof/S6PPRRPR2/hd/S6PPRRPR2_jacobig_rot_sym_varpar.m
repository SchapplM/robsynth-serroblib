% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRPR2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:33
	% EndTime: 2019-10-09 21:12:33
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t41 = sin(pkin(6));
	t44 = cos(pkin(7));
	t47 = t41 * t44;
	t42 = cos(pkin(12));
	t45 = cos(pkin(6));
	t46 = t42 * t45;
	t43 = cos(pkin(11));
	t40 = sin(pkin(7));
	t39 = sin(pkin(11));
	t38 = sin(pkin(12));
	t1 = [0, 0, -(-t43 * t38 - t39 * t46) * t40 + t39 * t47, 0, 0, 0; 0, 0, -(-t39 * t38 + t43 * t46) * t40 - t43 * t47, 0, 0, 0; 0, 0, -t41 * t42 * t40 + t45 * t44, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->13), mult. (46->32), div. (0->0), fcn. (67->10), ass. (0->17)
	t89 = sin(pkin(11));
	t95 = cos(pkin(6));
	t101 = t89 * t95;
	t90 = sin(pkin(7));
	t91 = sin(pkin(6));
	t100 = t90 * t91;
	t94 = cos(pkin(7));
	t99 = t91 * t94;
	t93 = cos(pkin(11));
	t98 = t93 * t95;
	t97 = cos(qJ(3));
	t96 = sin(qJ(3));
	t92 = cos(pkin(12));
	t88 = sin(pkin(12));
	t87 = -t92 * t101 - t93 * t88;
	t86 = -t89 * t88 + t92 * t98;
	t1 = [0, 0, -t87 * t90 + t89 * t99, (-t88 * t101 + t93 * t92) * t96 + (-t89 * t100 - t87 * t94) * t97, 0, 0; 0, 0, -t86 * t90 - t93 * t99, (t88 * t98 + t89 * t92) * t96 + (t93 * t100 - t86 * t94) * t97, 0, 0; 0, 0, -t92 * t100 + t95 * t94, -t95 * t90 * t97 + (-t92 * t94 * t97 + t88 * t96) * t91, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (15->13), mult. (46->32), div. (0->0), fcn. (67->10), ass. (0->17)
	t109 = sin(pkin(11));
	t115 = cos(pkin(6));
	t121 = t109 * t115;
	t110 = sin(pkin(7));
	t111 = sin(pkin(6));
	t120 = t110 * t111;
	t114 = cos(pkin(7));
	t119 = t111 * t114;
	t113 = cos(pkin(11));
	t118 = t113 * t115;
	t117 = cos(qJ(3));
	t116 = sin(qJ(3));
	t112 = cos(pkin(12));
	t108 = sin(pkin(12));
	t107 = -t113 * t108 - t112 * t121;
	t106 = -t109 * t108 + t112 * t118;
	t1 = [0, 0, -t107 * t110 + t109 * t119, (-t108 * t121 + t113 * t112) * t116 + (-t107 * t114 - t109 * t120) * t117, 0, 0; 0, 0, -t106 * t110 - t113 * t119, (t108 * t118 + t109 * t112) * t116 + (-t106 * t114 + t113 * t120) * t117, 0, 0; 0, 0, -t112 * t120 + t115 * t114, -t115 * t110 * t117 + (-t112 * t114 * t117 + t108 * t116) * t111, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:12:34
	% EndTime: 2019-10-09 21:12:34
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->20), mult. (98->46), div. (0->0), fcn. (140->12), ass. (0->28)
	t141 = sin(pkin(11));
	t147 = cos(pkin(6));
	t159 = t141 * t147;
	t142 = sin(pkin(7));
	t143 = sin(pkin(6));
	t158 = t142 * t143;
	t157 = t142 * t147;
	t146 = cos(pkin(7));
	t156 = t143 * t146;
	t144 = cos(pkin(12));
	t155 = t144 * t146;
	t145 = cos(pkin(11));
	t154 = t145 * t147;
	t140 = sin(pkin(12));
	t136 = -t141 * t140 + t144 * t154;
	t153 = -t136 * t146 + t145 * t158;
	t138 = -t145 * t140 - t144 * t159;
	t152 = t138 * t146 + t141 * t158;
	t151 = cos(qJ(3));
	t150 = cos(qJ(4));
	t149 = sin(qJ(3));
	t148 = sin(qJ(4));
	t139 = -t140 * t159 + t145 * t144;
	t137 = t140 * t154 + t141 * t144;
	t135 = -t144 * t158 + t147 * t146;
	t134 = -t138 * t142 + t141 * t156;
	t133 = -t136 * t142 - t145 * t156;
	t1 = [0, 0, t134, t139 * t149 - t152 * t151, 0, (t139 * t151 + t152 * t149) * t150 + t134 * t148; 0, 0, t133, t137 * t149 + t153 * t151, 0, (t137 * t151 - t153 * t149) * t150 + t133 * t148; 0, 0, t135, -t151 * t157 + (t140 * t149 - t151 * t155) * t143, 0, (t149 * t157 + (t140 * t151 + t149 * t155) * t143) * t150 + t135 * t148;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end