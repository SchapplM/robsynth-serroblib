% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:20
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRRRR2_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRRRR2_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t41 = sin(pkin(6));
	t44 = cos(pkin(7));
	t47 = t41 * t44;
	t42 = cos(pkin(13));
	t45 = cos(pkin(6));
	t46 = t42 * t45;
	t43 = cos(pkin(12));
	t40 = sin(pkin(7));
	t39 = sin(pkin(12));
	t38 = sin(pkin(13));
	t1 = [0, 0, -(-t43 * t38 - t39 * t46) * t40 + t39 * t47, 0, 0, 0; 0, 0, -(-t39 * t38 + t43 * t46) * t40 - t43 * t47, 0, 0, 0; 0, 0, -t41 * t42 * t40 + t45 * t44, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (15->13), mult. (46->32), div. (0->0), fcn. (67->10), ass. (0->17)
	t89 = sin(pkin(12));
	t95 = cos(pkin(6));
	t101 = t89 * t95;
	t90 = sin(pkin(7));
	t91 = sin(pkin(6));
	t100 = t90 * t91;
	t94 = cos(pkin(7));
	t99 = t91 * t94;
	t93 = cos(pkin(12));
	t98 = t93 * t95;
	t97 = cos(qJ(3));
	t96 = sin(qJ(3));
	t92 = cos(pkin(13));
	t88 = sin(pkin(13));
	t87 = -t92 * t101 - t93 * t88;
	t86 = -t89 * t88 + t92 * t98;
	t1 = [0, 0, -t87 * t90 + t89 * t99, (-t88 * t101 + t93 * t92) * t96 + (-t89 * t100 - t87 * t94) * t97, 0, 0; 0, 0, -t86 * t90 - t93 * t99, (t88 * t98 + t89 * t92) * t96 + (t93 * t100 - t86 * t94) * t97, 0, 0; 0, 0, -t92 * t100 + t95 * t94, -t95 * t90 * t97 + (-t92 * t94 * t97 + t88 * t96) * t91, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->20), mult. (98->46), div. (0->0), fcn. (140->12), ass. (0->28)
	t138 = sin(pkin(12));
	t144 = cos(pkin(6));
	t156 = t138 * t144;
	t139 = sin(pkin(7));
	t140 = sin(pkin(6));
	t155 = t139 * t140;
	t154 = t139 * t144;
	t143 = cos(pkin(7));
	t153 = t140 * t143;
	t141 = cos(pkin(13));
	t152 = t141 * t143;
	t142 = cos(pkin(12));
	t151 = t142 * t144;
	t137 = sin(pkin(13));
	t133 = -t138 * t137 + t141 * t151;
	t150 = -t133 * t143 + t142 * t155;
	t135 = -t142 * t137 - t141 * t156;
	t149 = t135 * t143 + t138 * t155;
	t148 = cos(qJ(3));
	t147 = cos(qJ(4));
	t146 = sin(qJ(3));
	t145 = sin(qJ(4));
	t136 = -t137 * t156 + t142 * t141;
	t134 = t137 * t151 + t138 * t141;
	t132 = -t141 * t155 + t144 * t143;
	t131 = -t135 * t139 + t138 * t153;
	t130 = -t133 * t139 - t142 * t153;
	t1 = [0, 0, t131, t136 * t146 - t149 * t148, (t136 * t148 + t149 * t146) * t145 - t131 * t147, 0; 0, 0, t130, t134 * t146 + t150 * t148, (t134 * t148 - t150 * t146) * t145 - t130 * t147, 0; 0, 0, t132, -t148 * t154 + (t137 * t146 - t148 * t152) * t140, (t146 * t154 + (t137 * t148 + t146 * t152) * t140) * t145 - t132 * t147, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:20:19
	% EndTime: 2019-10-09 21:20:19
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (51->20), mult. (150->46), div. (0->0), fcn. (213->12), ass. (0->31)
	t157 = sin(pkin(12));
	t163 = cos(pkin(6));
	t175 = t157 * t163;
	t158 = sin(pkin(7));
	t159 = sin(pkin(6));
	t174 = t158 * t159;
	t173 = t158 * t163;
	t162 = cos(pkin(7));
	t172 = t159 * t162;
	t160 = cos(pkin(13));
	t171 = t160 * t162;
	t161 = cos(pkin(12));
	t170 = t161 * t163;
	t156 = sin(pkin(13));
	t152 = -t157 * t156 + t160 * t170;
	t169 = -t152 * t162 + t161 * t174;
	t154 = -t161 * t156 - t160 * t175;
	t168 = t154 * t162 + t157 * t174;
	t167 = cos(qJ(3));
	t166 = cos(qJ(4));
	t165 = sin(qJ(3));
	t164 = sin(qJ(4));
	t155 = -t156 * t175 + t161 * t160;
	t153 = t156 * t170 + t157 * t160;
	t151 = -t160 * t174 + t163 * t162;
	t150 = -t154 * t158 + t157 * t172;
	t149 = -t152 * t158 - t161 * t172;
	t148 = (t165 * t173 + (t156 * t167 + t165 * t171) * t159) * t164 - t151 * t166;
	t147 = (t155 * t167 + t168 * t165) * t164 - t150 * t166;
	t146 = (t153 * t167 - t169 * t165) * t164 - t149 * t166;
	t1 = [0, 0, t150, t155 * t165 - t168 * t167, t147, t147; 0, 0, t149, t153 * t165 + t169 * t167, t146, t146; 0, 0, t151, -t167 * t173 + (t156 * t165 - t167 * t171) * t159, t148, t148;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end