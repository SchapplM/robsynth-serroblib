% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:08
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PPRPRR1_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PPRPRR1_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
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
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:43
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t56 = sin(pkin(6));
	t59 = cos(pkin(7));
	t62 = t56 * t59;
	t57 = cos(pkin(12));
	t60 = cos(pkin(6));
	t61 = t57 * t60;
	t58 = cos(pkin(11));
	t55 = sin(pkin(7));
	t54 = sin(pkin(11));
	t53 = sin(pkin(12));
	t1 = [0, 0, -(-t58 * t53 - t54 * t61) * t55 + t54 * t62, 0, 0, 0; 0, 0, -(-t54 * t53 + t58 * t61) * t55 - t58 * t62, 0, 0, 0; 0, 0, -t56 * t57 * t55 + t60 * t59, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:43
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (24->15), mult. (70->35), div. (0->0), fcn. (100->12), ass. (0->23)
	t108 = sin(pkin(11));
	t115 = cos(pkin(6));
	t122 = t108 * t115;
	t109 = sin(pkin(7));
	t106 = sin(pkin(13));
	t111 = cos(pkin(13));
	t116 = sin(qJ(3));
	t117 = cos(qJ(3));
	t118 = -t106 * t116 + t111 * t117;
	t101 = t118 * t109;
	t110 = sin(pkin(6));
	t121 = t110 * t101;
	t114 = cos(pkin(7));
	t120 = t110 * t114;
	t113 = cos(pkin(11));
	t119 = t113 * t115;
	t112 = cos(pkin(12));
	t107 = sin(pkin(12));
	t105 = -t117 * t106 - t116 * t111;
	t104 = -t113 * t107 - t112 * t122;
	t103 = -t108 * t107 + t112 * t119;
	t102 = t118 * t114;
	t1 = [0, 0, -t104 * t109 + t108 * t120, 0, -(-t107 * t122 + t113 * t112) * t105 - t104 * t102 - t108 * t121, 0; 0, 0, -t103 * t109 - t113 * t120, 0, -(t107 * t119 + t108 * t112) * t105 - t103 * t102 + t113 * t121, 0; 0, 0, -t110 * t112 * t109 + t115 * t114, 0, -t115 * t101 + (-t102 * t112 - t105 * t107) * t110, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:08:44
	% EndTime: 2019-10-09 21:08:44
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (51->24), mult. (146->54), div. (0->0), fcn. (206->14), ass. (0->33)
	t153 = sin(pkin(11));
	t155 = sin(pkin(6));
	t170 = t153 * t155;
	t160 = cos(pkin(6));
	t169 = t153 * t160;
	t158 = cos(pkin(11));
	t168 = t155 * t158;
	t159 = cos(pkin(7));
	t167 = t155 * t159;
	t166 = t158 * t160;
	t151 = sin(pkin(13));
	t156 = cos(pkin(13));
	t162 = sin(qJ(3));
	t164 = cos(qJ(3));
	t165 = t164 * t151 + t162 * t156;
	t150 = -t162 * t151 + t164 * t156;
	t163 = cos(qJ(5));
	t161 = sin(qJ(5));
	t157 = cos(pkin(12));
	t154 = sin(pkin(7));
	t152 = sin(pkin(12));
	t148 = -t152 * t169 + t158 * t157;
	t147 = -t158 * t152 - t157 * t169;
	t146 = t152 * t166 + t153 * t157;
	t145 = -t153 * t152 + t157 * t166;
	t144 = -t155 * t157 * t154 + t160 * t159;
	t143 = t165 * t159;
	t142 = t150 * t159;
	t141 = t165 * t154;
	t140 = t150 * t154;
	t139 = -t147 * t154 + t153 * t167;
	t138 = -t145 * t154 - t158 * t167;
	t1 = [0, 0, t139, 0, -t140 * t170 - t147 * t142 + t148 * t165, (t141 * t170 + t147 * t143 + t148 * t150) * t161 - t139 * t163; 0, 0, t138, 0, t140 * t168 - t145 * t142 + t146 * t165, (-t141 * t168 + t145 * t143 + t146 * t150) * t161 - t138 * t163; 0, 0, t144, 0, -t160 * t140 + (-t142 * t157 + t152 * t165) * t155, (t160 * t141 + (t143 * t157 + t150 * t152) * t155) * t161 - t144 * t163;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end