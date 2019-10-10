% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:31
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPPR9_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR9_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR9_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:24
	% EndTime: 2019-10-10 11:31:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:25
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t103 = sin(qJ(2));
	t104 = sin(qJ(1));
	t111 = t103 * t104;
	t106 = cos(qJ(1));
	t110 = t103 * t106;
	t105 = cos(qJ(2));
	t109 = t104 * t105;
	t108 = t105 * t106;
	t101 = sin(pkin(6));
	t107 = qJD(1) * t101;
	t102 = cos(pkin(6));
	t1 = [0, t106 * t107, (-t102 * t111 + t108) * qJD(2) + (t102 * t108 - t111) * qJD(1), 0, 0, 0; 0, t104 * t107, (t102 * t110 + t109) * qJD(2) + (t102 * t109 + t110) * qJD(1), 0, 0, 0; 0, 0, t101 * qJD(2) * t103, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:25
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t133 = sin(qJ(2));
	t134 = sin(qJ(1));
	t141 = t133 * t134;
	t136 = cos(qJ(1));
	t140 = t133 * t136;
	t135 = cos(qJ(2));
	t139 = t134 * t135;
	t138 = t135 * t136;
	t131 = sin(pkin(6));
	t137 = qJD(1) * t131;
	t132 = cos(pkin(6));
	t1 = [0, t136 * t137, (-t132 * t141 + t138) * qJD(2) + (t132 * t138 - t141) * qJD(1), 0, 0, 0; 0, t134 * t137, (t132 * t140 + t139) * qJD(2) + (t132 * t139 + t140) * qJD(1), 0, 0, 0; 0, 0, t131 * qJD(2) * t133, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:25
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t158 = sin(qJ(2));
	t159 = sin(qJ(1));
	t166 = t158 * t159;
	t161 = cos(qJ(1));
	t165 = t158 * t161;
	t160 = cos(qJ(2));
	t164 = t159 * t160;
	t163 = t160 * t161;
	t156 = sin(pkin(6));
	t162 = qJD(1) * t156;
	t157 = cos(pkin(6));
	t1 = [0, t161 * t162, (-t157 * t166 + t163) * qJD(2) + (t157 * t163 - t166) * qJD(1), 0, 0, 0; 0, t159 * t162, (t157 * t165 + t164) * qJD(2) + (t157 * t164 + t165) * qJD(1), 0, 0, 0; 0, 0, t156 * qJD(2) * t158, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:31:25
	% EndTime: 2019-10-10 11:31:25
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
	t203 = sin(pkin(6));
	t208 = cos(qJ(3));
	t222 = t203 * t208;
	t210 = cos(qJ(1));
	t221 = t203 * t210;
	t206 = sin(qJ(2));
	t207 = sin(qJ(1));
	t220 = t206 * t207;
	t219 = t206 * t210;
	t209 = cos(qJ(2));
	t218 = t207 * t209;
	t217 = t210 * t209;
	t216 = qJD(1) * t203;
	t205 = sin(qJ(3));
	t215 = qJD(2) * t205;
	t204 = cos(pkin(6));
	t214 = t204 * t217 - t220;
	t213 = t204 * t218 + t219;
	t212 = t204 * t219 + t218;
	t211 = t204 * t220 - t217;
	t1 = [0, t210 * t216, t214 * qJD(1) - t211 * qJD(2), 0, 0, (-t207 * t203 * t205 + t211 * t208) * qJD(3) + t213 * t215 + (t212 * t205 + t208 * t221) * qJD(1); 0, t207 * t216, t213 * qJD(1) + t212 * qJD(2), 0, 0, (t205 * t221 - t212 * t208) * qJD(3) - t214 * t215 + (t211 * t205 + t207 * t222) * qJD(1); 0, 0, t203 * qJD(2) * t206, 0, 0, -t203 * t209 * t215 + (-t204 * t205 - t206 * t222) * qJD(3);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end