% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPRRR2
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
% Datum: 2019-10-09 21:55
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPRRR2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:55:41
	% EndTime: 2019-10-09 21:55:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:55:41
	% EndTime: 2019-10-09 21:55:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:55:41
	% EndTime: 2019-10-09 21:55:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:55:41
	% EndTime: 2019-10-09 21:55:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:55:41
	% EndTime: 2019-10-09 21:55:42
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
	% StartTime: 2019-10-09 21:55:42
	% EndTime: 2019-10-09 21:55:42
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (27->14), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->19)
	t172 = sin(pkin(12));
	t175 = cos(pkin(12));
	t179 = sin(qJ(2));
	t181 = cos(qJ(2));
	t183 = t179 * t172 - t181 * t175;
	t186 = t183 * qJD(2);
	t184 = t172 * t181 + t175 * t179;
	t170 = t184 * qJD(2);
	t174 = sin(pkin(6));
	t178 = sin(qJ(4));
	t185 = t174 * t178;
	t180 = cos(qJ(4));
	t177 = cos(pkin(6));
	t176 = cos(pkin(11));
	t173 = sin(pkin(11));
	t168 = t184 * t177;
	t167 = t177 * t186;
	t166 = t177 * t170;
	t1 = [0, 0, 0, -t173 * t166 - t176 * t186, (t173 * t167 - t176 * t170) * t178 + ((-t173 * t168 - t176 * t183) * t180 + t173 * t185) * qJD(4), 0; 0, 0, 0, t176 * t166 - t173 * t186, (-t176 * t167 - t173 * t170) * t178 + ((t176 * t168 - t173 * t183) * t180 - t176 * t185) * qJD(4), 0; 0, 0, 0, t174 * t170, t177 * qJD(4) * t178 + (t184 * qJD(4) * t180 - t178 * t186) * t174, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:55:42
	% EndTime: 2019-10-09 21:55:42
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (47->14), mult. (164->37), div. (0->0), fcn. (178->10), ass. (0->22)
	t201 = sin(pkin(12));
	t204 = cos(pkin(12));
	t208 = sin(qJ(2));
	t210 = cos(qJ(2));
	t212 = t208 * t201 - t210 * t204;
	t215 = t212 * qJD(2);
	t213 = t201 * t210 + t204 * t208;
	t199 = t213 * qJD(2);
	t203 = sin(pkin(6));
	t207 = sin(qJ(4));
	t214 = t203 * t207;
	t209 = cos(qJ(4));
	t206 = cos(pkin(6));
	t205 = cos(pkin(11));
	t202 = sin(pkin(11));
	t197 = t213 * t206;
	t196 = t206 * t215;
	t195 = t206 * t199;
	t194 = t206 * qJD(4) * t207 + (t213 * qJD(4) * t209 - t207 * t215) * t203;
	t193 = (t202 * t196 - t205 * t199) * t207 + ((-t202 * t197 - t205 * t212) * t209 + t202 * t214) * qJD(4);
	t192 = (-t205 * t196 - t202 * t199) * t207 + ((t205 * t197 - t202 * t212) * t209 - t205 * t214) * qJD(4);
	t1 = [0, 0, 0, -t202 * t195 - t205 * t215, t193, t193; 0, 0, 0, t205 * t195 - t202 * t215, t192, t192; 0, 0, 0, t203 * t199, t194, t194;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end