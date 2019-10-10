% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:10
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRPR1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRPR1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
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
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->13), mult. (46->32), div. (0->0), fcn. (67->10), ass. (0->17)
	t123 = sin(pkin(11));
	t129 = cos(pkin(6));
	t135 = t123 * t129;
	t124 = sin(pkin(7));
	t125 = sin(pkin(6));
	t134 = t124 * t125;
	t128 = cos(pkin(7));
	t133 = t125 * t128;
	t127 = cos(pkin(11));
	t132 = t127 * t129;
	t131 = cos(qJ(3));
	t130 = sin(qJ(3));
	t126 = cos(pkin(12));
	t122 = sin(pkin(12));
	t121 = -t127 * t122 - t126 * t135;
	t120 = -t123 * t122 + t126 * t132;
	t1 = [0, 0, -t121 * t124 + t123 * t133, (-t122 * t135 + t127 * t126) * t130 + (-t121 * t128 - t123 * t134) * t131, 0, 0; 0, 0, -t120 * t124 - t127 * t133, (t122 * t132 + t123 * t126) * t130 + (-t120 * t128 + t127 * t134) * t131, 0, 0; 0, 0, -t126 * t134 + t129 * t128, -t129 * t124 * t131 + (-t126 * t128 * t131 + t122 * t130) * t125, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:10:40
	% EndTime: 2019-10-09 21:10:40
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (33->20), mult. (98->46), div. (0->0), fcn. (140->12), ass. (0->28)
	t148 = sin(pkin(11));
	t154 = cos(pkin(6));
	t166 = t148 * t154;
	t149 = sin(pkin(7));
	t150 = sin(pkin(6));
	t165 = t149 * t150;
	t164 = t149 * t154;
	t153 = cos(pkin(7));
	t163 = t150 * t153;
	t151 = cos(pkin(12));
	t162 = t151 * t153;
	t152 = cos(pkin(11));
	t161 = t152 * t154;
	t147 = sin(pkin(12));
	t143 = -t148 * t147 + t151 * t161;
	t160 = -t143 * t153 + t152 * t165;
	t145 = -t152 * t147 - t151 * t166;
	t159 = t145 * t153 + t148 * t165;
	t158 = cos(qJ(3));
	t157 = cos(qJ(4));
	t156 = sin(qJ(3));
	t155 = sin(qJ(4));
	t146 = -t147 * t166 + t152 * t151;
	t144 = t147 * t161 + t148 * t151;
	t142 = -t151 * t165 + t154 * t153;
	t141 = -t145 * t149 + t148 * t163;
	t140 = -t143 * t149 - t152 * t163;
	t1 = [0, 0, t141, t146 * t156 - t159 * t158, 0, (t146 * t158 + t159 * t156) * t155 - t141 * t157; 0, 0, t140, t144 * t156 + t160 * t158, 0, (t144 * t158 - t160 * t156) * t155 - t140 * t157; 0, 0, t142, -t158 * t164 + (t147 * t156 - t158 * t162) * t150, 0, (t156 * t164 + (t147 * t158 + t156 * t162) * t150) * t155 - t142 * t157;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end