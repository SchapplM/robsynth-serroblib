% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 21:26
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPPRR2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRPPRR2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (7->4), mult. (30->12), div. (0->0), fcn. (30->8), ass. (0->10)
	t109 = sin(pkin(11));
	t112 = cos(pkin(11));
	t115 = sin(qJ(2));
	t116 = cos(qJ(2));
	t117 = qJD(2) * (-t109 * t115 + t112 * t116);
	t113 = cos(pkin(10));
	t110 = sin(pkin(10));
	t108 = (-t109 * t116 - t112 * t115) * qJD(2);
	t107 = cos(pkin(6)) * t117;
	t1 = [0, 0, 0, 0, -t110 * t107 + t113 * t108, 0; 0, 0, 0, 0, t113 * t107 + t110 * t108, 0; 0, 0, 0, 0, sin(pkin(6)) * t117, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 21:26:19
	% EndTime: 2019-10-09 21:26:19
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (27->15), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->19)
	t175 = sin(pkin(11));
	t178 = cos(pkin(11));
	t182 = sin(qJ(2));
	t184 = cos(qJ(2));
	t186 = t175 * t182 - t178 * t184;
	t189 = t186 * qJD(2);
	t187 = t184 * t175 + t182 * t178;
	t173 = t187 * qJD(2);
	t177 = sin(pkin(6));
	t183 = cos(qJ(5));
	t188 = t177 * t183;
	t181 = sin(qJ(5));
	t180 = cos(pkin(6));
	t179 = cos(pkin(10));
	t176 = sin(pkin(10));
	t171 = t186 * t180;
	t170 = t180 * t189;
	t169 = t180 * t173;
	t1 = [0, 0, 0, 0, t176 * t170 - t179 * t173, -(-t176 * t169 - t179 * t189) * t183 + (t176 * t188 + (-t176 * t171 + t179 * t187) * t181) * qJD(5); 0, 0, 0, 0, -t179 * t170 - t176 * t173, -(t179 * t169 - t176 * t189) * t183 + (-t179 * t188 + (t179 * t171 + t176 * t187) * t181) * qJD(5); 0, 0, 0, 0, -t177 * t189, t180 * qJD(5) * t183 + (t186 * qJD(5) * t181 - t183 * t173) * t177;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end