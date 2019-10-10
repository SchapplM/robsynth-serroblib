% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RRPR2
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:51
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RRPR2_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RRPR2_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR2_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_jacobiaD_transl_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:51
	% EndTime: 2019-10-09 20:51:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:51
	% EndTime: 2019-10-09 20:51:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (8->6), div. (0->0), fcn. (4->2), ass. (0->3)
	t27 = cos(qJ(1));
	t26 = sin(qJ(1));
	t1 = [(-r_i_i_C(1) * t27 + r_i_i_C(2) * t26) * qJD(1), 0, 0, 0; (-r_i_i_C(1) * t26 - r_i_i_C(2) * t27) * qJD(1), 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:51
	% EndTime: 2019-10-09 20:51:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (22->6), mult. (20->9), div. (0->0), fcn. (10->4), ass. (0->8)
	t43 = pkin(1) * qJD(1);
	t40 = qJ(1) + qJ(2);
	t37 = sin(t40);
	t38 = cos(t40);
	t39 = qJD(1) + qJD(2);
	t42 = (-r_i_i_C(1) * t38 + r_i_i_C(2) * t37) * t39;
	t41 = (-r_i_i_C(1) * t37 - r_i_i_C(2) * t38) * t39;
	t1 = [-cos(qJ(1)) * t43 + t42, t42, 0, 0; -sin(qJ(1)) * t43 + t41, t41, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:51
	% EndTime: 2019-10-09 20:51:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (58->10), mult. (42->12), div. (0->0), fcn. (24->4), ass. (0->12)
	t31 = qJ(3) + r_i_i_C(3);
	t30 = -pkin(2) - r_i_i_C(1);
	t24 = qJ(1) + qJ(2);
	t21 = sin(t24);
	t23 = qJD(1) + qJD(2);
	t29 = t23 * t21;
	t22 = cos(t24);
	t28 = t23 * t22;
	t27 = pkin(1) * qJD(1);
	t26 = t21 * qJD(3) + t31 * t28 + t30 * t29;
	t25 = t22 * qJD(3) + (-t31 * t21 + t30 * t22) * t23;
	t1 = [-cos(qJ(1)) * t27 + t25, t25, t28, 0; -sin(qJ(1)) * t27 + t26, t26, t29, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:51
	% EndTime: 2019-10-09 20:51:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (164->19), mult. (142->27), div. (0->0), fcn. (116->6), ass. (0->19)
	t103 = qJ(1) + qJ(2);
	t100 = sin(t103);
	t101 = cos(t103);
	t102 = qJD(1) + qJD(2);
	t104 = sin(qJ(4));
	t105 = cos(qJ(4));
	t108 = qJD(4) * t105;
	t109 = qJD(4) * t104;
	t91 = -t100 * t109 - t101 * t108 + (t100 * t104 + t101 * t105) * t102;
	t110 = t102 * t101;
	t111 = t102 * t100;
	t92 = t100 * t108 - t101 * t109 + t104 * t110 - t105 * t111;
	t115 = t91 * r_i_i_C(1) - t92 * r_i_i_C(2);
	t114 = -t92 * r_i_i_C(1) - t91 * r_i_i_C(2);
	t113 = -pkin(2) - pkin(3);
	t112 = pkin(1) * qJD(1);
	t107 = qJ(3) * t110 + t100 * qJD(3) + t113 * t111 - t114;
	t106 = t101 * qJD(3) + (-qJ(3) * t100 + t113 * t101) * t102 - t115;
	t1 = [-cos(qJ(1)) * t112 + t106, t106, t110, t115; -sin(qJ(1)) * t112 + t107, t107, t111, t114; 0, 0, 0, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end