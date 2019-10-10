% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRR5
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:22
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRR5_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR5_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRR5_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (12->6), mult. (48->16), div. (0->0), fcn. (48->6), ass. (0->15)
	t128 = sin(qJ(2));
	t129 = sin(qJ(1));
	t136 = t128 * t129;
	t131 = cos(qJ(1));
	t135 = t128 * t131;
	t130 = cos(qJ(2));
	t134 = t129 * t130;
	t133 = t130 * t131;
	t126 = sin(pkin(6));
	t132 = qJD(1) * t126;
	t127 = cos(pkin(6));
	t125 = t126 * qJD(2) * t128;
	t124 = (t127 * t135 + t134) * qJD(2) + (t127 * t134 + t135) * qJD(1);
	t123 = (-t127 * t136 + t133) * qJD(2) + (t127 * t133 - t136) * qJD(1);
	t1 = [0, t131 * t132, t123, t123, 0, 0; 0, t129 * t132, t124, t124, 0, 0; 0, 0, t125, t125, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:26
	% EndTime: 2019-10-10 13:22:26
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (18->6), mult. (70->16), div. (0->0), fcn. (70->6), ass. (0->15)
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
	t130 = t131 * qJD(2) * t133;
	t129 = (t132 * t140 + t139) * qJD(2) + (t132 * t139 + t140) * qJD(1);
	t128 = (-t132 * t141 + t138) * qJD(2) + (t132 * t138 - t141) * qJD(1);
	t1 = [0, t136 * t137, t128, t128, t128, 0; 0, t134 * t137, t129, t129, t129, 0; 0, 0, t130, t130, t130, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:22:27
	% EndTime: 2019-10-10 13:22:27
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (68->18), mult. (122->40), div. (0->0), fcn. (124->8), ass. (0->27)
	t202 = sin(pkin(6));
	t205 = sin(qJ(1));
	t220 = t202 * t205;
	t207 = cos(qJ(1));
	t219 = t202 * t207;
	t204 = sin(qJ(2));
	t218 = t204 * t205;
	t217 = t204 * t207;
	t206 = cos(qJ(2));
	t216 = t205 * t206;
	t215 = t207 * t206;
	t214 = qJD(1) * t202;
	t201 = qJ(3) + qJ(4) + qJ(5);
	t198 = sin(t201);
	t213 = qJD(2) * t198;
	t212 = qJD(2) * t202;
	t203 = cos(pkin(6));
	t211 = t203 * t215 - t218;
	t210 = t203 * t216 + t217;
	t209 = t203 * t217 + t216;
	t208 = -t203 * t218 + t215;
	t200 = qJD(3) + qJD(4) + qJD(5);
	t199 = cos(t201);
	t197 = t204 * t212;
	t196 = t210 * qJD(1) + t209 * qJD(2);
	t195 = t211 * qJD(1) + t208 * qJD(2);
	t1 = [0, t207 * t214, t195, t195, t195, (t198 * t220 + t208 * t199) * t200 - t210 * t213 + (-t209 * t198 - t199 * t219) * qJD(1); 0, t205 * t214, t196, t196, t196, (-t198 * t219 + t209 * t199) * t200 + t211 * t213 + (t208 * t198 - t199 * t220) * qJD(1); 0, 0, t197, t197, t197, t202 * t204 * t200 * t199 + (t200 * t203 + t206 * t212) * t198;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end