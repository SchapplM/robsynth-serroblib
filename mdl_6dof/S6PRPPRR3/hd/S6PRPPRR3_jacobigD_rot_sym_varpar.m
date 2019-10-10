% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:28
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPPRR3_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR3_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:08
	% EndTime: 2019-10-09 21:28:08
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (30->20), div. (0->0), fcn. (30->8), ass. (0->10)
	t114 = cos(pkin(6));
	t115 = sin(qJ(2));
	t118 = t114 * t115;
	t116 = cos(qJ(2));
	t117 = t114 * t116;
	t113 = cos(pkin(10));
	t112 = cos(pkin(11));
	t110 = sin(pkin(10));
	t109 = sin(pkin(11));
	t1 = [0, 0, 0, 0, ((-t110 * t117 - t113 * t115) * t109 - (-t110 * t118 + t113 * t116) * t112) * qJD(2), 0; 0, 0, 0, 0, ((-t110 * t115 + t113 * t117) * t109 - (t110 * t116 + t113 * t118) * t112) * qJD(2), 0; 0, 0, 0, 0, (t109 * t116 - t112 * t115) * sin(pkin(6)) * qJD(2), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:28:09
	% EndTime: 2019-10-09 21:28:09
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (27->18), mult. (97->48), div. (0->0), fcn. (104->10), ass. (0->23)
	t178 = sin(pkin(6));
	t182 = sin(qJ(5));
	t193 = t178 * t182;
	t181 = cos(pkin(6));
	t183 = sin(qJ(2));
	t192 = t181 * t183;
	t185 = cos(qJ(2));
	t191 = t181 * t185;
	t176 = sin(pkin(11));
	t179 = cos(pkin(11));
	t190 = t176 * t185 - t179 * t183;
	t177 = sin(pkin(10));
	t180 = cos(pkin(10));
	t189 = -t177 * t183 + t180 * t191;
	t188 = t177 * t185 + t180 * t192;
	t187 = t177 * t191 + t180 * t183;
	t186 = -t177 * t192 + t180 * t185;
	t184 = cos(qJ(5));
	t175 = t186 * qJD(2);
	t174 = t187 * qJD(2);
	t173 = t188 * qJD(2);
	t172 = t189 * qJD(2);
	t1 = [0, 0, 0, 0, -t174 * t176 - t175 * t179, (-t174 * t179 + t175 * t176) * t182 + ((t187 * t176 + t186 * t179) * t184 - t177 * t193) * qJD(5); 0, 0, 0, 0, t172 * t176 - t173 * t179, (t172 * t179 + t173 * t176) * t182 + ((-t189 * t176 + t188 * t179) * t184 + t180 * t193) * qJD(5); 0, 0, 0, 0, t190 * t178 * qJD(2), -t181 * qJD(5) * t182 + (-t190 * qJD(5) * t184 + (t183 * t176 + t185 * t179) * t182 * qJD(2)) * t178;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end