% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:48
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPR2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:33
	% EndTime: 2019-10-09 22:48:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = sin(qJ(2));
	t84 = cos(pkin(6)) * t82;
	t83 = cos(qJ(2));
	t80 = cos(pkin(11));
	t79 = sin(pkin(11));
	t1 = [0, 0, (-t79 * t84 + t80 * t83) * qJD(2), 0, 0, 0; 0, 0, (t79 * t83 + t80 * t84) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->2), mult. (24->9), div. (0->0), fcn. (24->6), ass. (0->9)
	t111 = sin(qJ(2));
	t113 = cos(pkin(6)) * t111;
	t112 = cos(qJ(2));
	t109 = cos(pkin(11));
	t108 = sin(pkin(11));
	t107 = sin(pkin(6)) * qJD(2) * t111;
	t106 = (-t108 * t113 + t109 * t112) * qJD(2);
	t105 = (t108 * t112 + t109 * t113) * qJD(2);
	t1 = [0, 0, t106, t106, 0, 0; 0, 0, t105, t105, 0, 0; 0, 0, t107, t107, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->2), mult. (24->9), div. (0->0), fcn. (24->6), ass. (0->9)
	t158 = sin(qJ(2));
	t160 = cos(pkin(6)) * t158;
	t159 = cos(qJ(2));
	t156 = cos(pkin(11));
	t155 = sin(pkin(11));
	t154 = sin(pkin(6)) * qJD(2) * t158;
	t153 = (-t155 * t160 + t156 * t159) * qJD(2);
	t152 = (t155 * t159 + t156 * t160) * qJD(2);
	t1 = [0, 0, t153, t153, 0, 0; 0, 0, t152, t152, 0, 0; 0, 0, t154, t154, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:48:34
	% EndTime: 2019-10-09 22:48:34
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (29->12), mult. (60->30), div. (0->0), fcn. (62->8), ass. (0->21)
	t179 = qJ(3) + qJ(4);
	t176 = sin(t179);
	t181 = sin(pkin(6));
	t192 = t181 * t176;
	t183 = cos(pkin(6));
	t184 = sin(qJ(2));
	t191 = t183 * t184;
	t185 = cos(qJ(2));
	t190 = t183 * t185;
	t189 = qJD(2) * t176;
	t188 = qJD(2) * t181;
	t180 = sin(pkin(11));
	t182 = cos(pkin(11));
	t187 = t180 * t185 + t182 * t191;
	t186 = -t180 * t191 + t182 * t185;
	t178 = qJD(3) + qJD(4);
	t177 = cos(t179);
	t175 = t184 * t188;
	t174 = t186 * qJD(2);
	t173 = t187 * qJD(2);
	t1 = [0, 0, t174, t174, 0, (t186 * t177 + t180 * t192) * t178 + (-t180 * t190 - t182 * t184) * t189; 0, 0, t173, t173, 0, (t187 * t177 - t182 * t192) * t178 + (-t180 * t184 + t182 * t190) * t189; 0, 0, t175, t175, 0, t181 * t184 * t178 * t177 + (t178 * t183 + t185 * t188) * t176;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end