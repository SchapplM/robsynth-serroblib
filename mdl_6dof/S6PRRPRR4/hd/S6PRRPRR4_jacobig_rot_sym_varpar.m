% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:31
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRR4_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR4_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobig_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t18 = sin(pkin(6));
	t1 = [0, sin(pkin(11)) * t18, 0, 0, 0, 0; 0, -cos(pkin(11)) * t18, 0, 0, 0, 0; 0, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t57 = cos(pkin(6));
	t59 = cos(qJ(2));
	t60 = t57 * t59;
	t58 = sin(qJ(2));
	t56 = cos(pkin(11));
	t55 = sin(pkin(6));
	t54 = sin(pkin(11));
	t1 = [0, t54 * t55, t54 * t60 + t56 * t58, 0, 0, 0; 0, -t56 * t55, t54 * t58 - t56 * t60, 0, 0, 0; 0, t57, -t55 * t59, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->4), mult. (9->8), div. (0->0), fcn. (17->6), ass. (0->8)
	t75 = cos(pkin(6));
	t77 = cos(qJ(2));
	t78 = t75 * t77;
	t76 = sin(qJ(2));
	t74 = cos(pkin(11));
	t73 = sin(pkin(6));
	t72 = sin(pkin(11));
	t1 = [0, t72 * t73, t72 * t78 + t74 * t76, 0, 0, 0; 0, -t74 * t73, t72 * t76 - t74 * t78, 0, 0, 0; 0, t75, -t73 * t77, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:29
	% EndTime: 2019-10-09 22:31:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (16->8), div. (0->0), fcn. (29->6), ass. (0->11)
	t74 = sin(pkin(6));
	t78 = cos(qJ(2));
	t80 = t74 * t78;
	t76 = cos(pkin(6));
	t79 = t76 * t78;
	t77 = sin(qJ(2));
	t75 = cos(pkin(11));
	t73 = sin(pkin(11));
	t72 = t73 * t79 + t75 * t77;
	t71 = t73 * t77 - t75 * t79;
	t1 = [0, t73 * t74, t72, 0, -t72, 0; 0, -t75 * t74, t71, 0, -t71, 0; 0, t76, -t80, 0, t80, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:31:30
	% EndTime: 2019-10-09 22:31:30
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (19->17), mult. (52->33), div. (0->0), fcn. (81->10), ass. (0->21)
	t152 = sin(pkin(6));
	t156 = sin(qJ(3));
	t166 = t152 * t156;
	t159 = cos(qJ(3));
	t165 = t152 * t159;
	t160 = cos(qJ(2));
	t164 = t152 * t160;
	t153 = cos(pkin(11));
	t163 = t153 * t152;
	t154 = cos(pkin(6));
	t157 = sin(qJ(2));
	t162 = t154 * t157;
	t161 = t154 * t160;
	t158 = cos(qJ(5));
	t155 = sin(qJ(5));
	t151 = sin(pkin(11));
	t150 = -t151 * t162 + t153 * t160;
	t149 = t151 * t161 + t153 * t157;
	t148 = t151 * t160 + t153 * t162;
	t147 = t151 * t157 - t153 * t161;
	t1 = [0, t151 * t152, t149, 0, -t149, (t150 * t159 + t151 * t166) * t155 - (t150 * t156 - t151 * t165) * t158; 0, -t163, t147, 0, -t147, (t148 * t159 - t156 * t163) * t155 - (t148 * t156 + t159 * t163) * t158; 0, t154, -t164, 0, t164, (t154 * t156 + t157 * t165) * t155 - (-t154 * t159 + t157 * t166) * t158;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end