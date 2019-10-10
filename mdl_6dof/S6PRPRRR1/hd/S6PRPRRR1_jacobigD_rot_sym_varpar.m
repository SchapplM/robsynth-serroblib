% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:53
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:47
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (30->12), div. (0->0), fcn. (30->8), ass. (0->10)
	t107 = sin(pkin(12));
	t110 = cos(pkin(12));
	t113 = sin(qJ(2));
	t114 = cos(qJ(2));
	t116 = qJD(2) * (t107 * t114 + t110 * t113);
	t111 = cos(pkin(11));
	t108 = sin(pkin(11));
	t106 = (t107 * t113 - t110 * t114) * qJD(2);
	t105 = cos(pkin(6)) * t116;
	t1 = [0, 0, 0, -t108 * t105 - t111 * t106, 0, 0; 0, 0, 0, t111 * t105 - t108 * t106, 0, 0; 0, 0, 0, sin(pkin(6)) * t116, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:48
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (14->4), mult. (60->12), div. (0->0), fcn. (60->8), ass. (0->13)
	t135 = sin(pkin(12));
	t138 = cos(pkin(12));
	t141 = sin(qJ(2));
	t142 = cos(qJ(2));
	t144 = qJD(2) * (t135 * t142 + t138 * t141);
	t139 = cos(pkin(11));
	t136 = sin(pkin(11));
	t134 = (t135 * t141 - t138 * t142) * qJD(2);
	t133 = cos(pkin(6)) * t144;
	t132 = sin(pkin(6)) * t144;
	t131 = -t136 * t133 - t139 * t134;
	t130 = t139 * t133 - t136 * t134;
	t1 = [0, 0, 0, t131, t131, 0; 0, 0, 0, t130, t130, 0; 0, 0, 0, t132, t132, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:53:48
	% EndTime: 2019-10-09 21:53:49
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (49->15), mult. (127->33), div. (0->0), fcn. (134->10), ass. (0->23)
	t216 = sin(pkin(12));
	t219 = cos(pkin(12));
	t222 = sin(qJ(2));
	t223 = cos(qJ(2));
	t225 = t222 * t216 - t223 * t219;
	t229 = t225 * qJD(2);
	t226 = t216 * t223 + t219 * t222;
	t210 = t226 * qJD(2);
	t214 = qJD(4) + qJD(5);
	t215 = qJ(4) + qJ(5);
	t228 = cos(t215) * t214;
	t218 = sin(pkin(6));
	t221 = cos(pkin(6));
	t227 = t214 * t218 + t221 * t229;
	t220 = cos(pkin(11));
	t217 = sin(pkin(11));
	t212 = sin(t215);
	t208 = t226 * t221;
	t206 = t221 * t210;
	t205 = t218 * t210;
	t204 = -t217 * t206 - t220 * t229;
	t203 = t220 * t206 - t217 * t229;
	t1 = [0, 0, 0, t204, t204, (-t217 * t208 - t220 * t225) * t228 + (-t220 * t210 + t227 * t217) * t212; 0, 0, 0, t203, t203, (t220 * t208 - t217 * t225) * t228 + (-t217 * t210 - t227 * t220) * t212; 0, 0, 0, t205, t205, t221 * t214 * t212 + (-t212 * t229 + t226 * t228) * t218;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end