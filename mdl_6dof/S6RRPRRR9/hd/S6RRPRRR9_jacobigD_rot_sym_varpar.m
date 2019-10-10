% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:03
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRR9_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR9_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRR9_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t81 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t81, 0, 0, 0, 0; 0, sin(qJ(1)) * t81, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t108 = sin(qJ(2));
	t109 = sin(qJ(1));
	t116 = t108 * t109;
	t111 = cos(qJ(1));
	t115 = t108 * t111;
	t110 = cos(qJ(2));
	t114 = t109 * t110;
	t113 = t110 * t111;
	t106 = sin(pkin(6));
	t112 = qJD(1) * t106;
	t107 = cos(pkin(6));
	t1 = [0, t111 * t112, 0, (-t107 * t116 + t113) * qJD(2) + (t107 * t113 - t116) * qJD(1), 0, 0; 0, t109 * t112, 0, (t107 * t115 + t114) * qJD(2) + (t107 * t114 + t115) * qJD(1), 0, 0; 0, 0, 0, t106 * qJD(2) * t108, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:13
	% EndTime: 2019-10-10 11:03:13
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (12->6), mult. (48->16), div. (0->0), fcn. (48->6), ass. (0->15)
	t130 = sin(qJ(2));
	t131 = sin(qJ(1));
	t138 = t130 * t131;
	t133 = cos(qJ(1));
	t137 = t130 * t133;
	t132 = cos(qJ(2));
	t136 = t131 * t132;
	t135 = t132 * t133;
	t128 = sin(pkin(6));
	t134 = qJD(1) * t128;
	t129 = cos(pkin(6));
	t127 = t128 * qJD(2) * t130;
	t126 = (t129 * t137 + t136) * qJD(2) + (t129 * t136 + t137) * qJD(1);
	t125 = (-t129 * t138 + t135) * qJD(2) + (t129 * t135 - t138) * qJD(1);
	t1 = [0, t133 * t134, 0, t125, t125, 0; 0, t131 * t134, 0, t126, t126, 0; 0, 0, 0, t127, t127, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:03:14
	% EndTime: 2019-10-10 11:03:14
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (56->18), mult. (100->40), div. (0->0), fcn. (102->8), ass. (0->27)
	t199 = sin(pkin(6));
	t202 = sin(qJ(1));
	t217 = t199 * t202;
	t204 = cos(qJ(1));
	t216 = t199 * t204;
	t201 = sin(qJ(2));
	t215 = t201 * t202;
	t214 = t201 * t204;
	t203 = cos(qJ(2));
	t213 = t202 * t203;
	t212 = t204 * t203;
	t211 = qJD(1) * t199;
	t197 = pkin(12) + qJ(4) + qJ(5);
	t195 = sin(t197);
	t210 = qJD(2) * t195;
	t209 = qJD(2) * t199;
	t200 = cos(pkin(6));
	t208 = t200 * t212 - t215;
	t207 = t200 * t213 + t214;
	t206 = t200 * t214 + t213;
	t205 = -t200 * t215 + t212;
	t198 = qJD(4) + qJD(5);
	t196 = cos(t197);
	t194 = t201 * t209;
	t193 = qJD(1) * t207 + qJD(2) * t206;
	t192 = qJD(1) * t208 + qJD(2) * t205;
	t1 = [0, t204 * t211, 0, t192, t192, (t195 * t217 + t196 * t205) * t198 - t207 * t210 + (-t195 * t206 - t196 * t216) * qJD(1); 0, t202 * t211, 0, t193, t193, (-t195 * t216 + t196 * t206) * t198 + t208 * t210 + (t195 * t205 - t196 * t217) * qJD(1); 0, 0, 0, t194, t194, t199 * t201 * t198 * t196 + (t198 * t200 + t203 * t209) * t195;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end