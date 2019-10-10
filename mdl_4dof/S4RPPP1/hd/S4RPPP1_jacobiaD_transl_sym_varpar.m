% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S4RPPP1
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
% 
% Output:
% JaD_transl [3x4]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:39
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S4RPPP1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_jacobiaD_transl_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_jacobiaD_transl_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S4RPPP1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPP1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_jacobiaD_transl_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:01
	% EndTime: 2019-10-09 20:39:01
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
	% StartTime: 2019-10-09 20:39:02
	% EndTime: 2019-10-09 20:39:02
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (14->13), mult. (48->27), div. (0->0), fcn. (38->6), ass. (0->12)
	t114 = sin(pkin(4));
	t124 = t114 * (r_i_i_C(3) + qJ(2));
	t116 = cos(pkin(4));
	t117 = sin(qJ(1));
	t122 = t116 * t117;
	t118 = cos(qJ(1));
	t121 = t116 * t118;
	t120 = qJD(1) * t114;
	t119 = t114 * qJD(2);
	t115 = cos(pkin(6));
	t113 = sin(pkin(6));
	t1 = [t118 * t119 + ((t113 * t122 - t115 * t118) * r_i_i_C(1) + (t113 * t118 + t115 * t122) * r_i_i_C(2) - t118 * pkin(1) - t117 * t124) * qJD(1), t118 * t120, 0, 0; t117 * t119 + ((-t113 * t121 - t115 * t117) * r_i_i_C(1) + (t113 * t117 - t115 * t121) * r_i_i_C(2) - t117 * pkin(1) + t118 * t124) * qJD(1), t117 * t120, 0, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:02
	% EndTime: 2019-10-09 20:39:02
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (28->19), mult. (90->32), div. (0->0), fcn. (78->6), ass. (0->21)
	t158 = r_i_i_C(1) + qJ(2);
	t157 = -r_i_i_C(3) - qJ(3);
	t142 = sin(pkin(4));
	t146 = cos(qJ(1));
	t156 = t142 * t146;
	t143 = cos(pkin(6));
	t145 = sin(qJ(1));
	t155 = t143 * t145;
	t154 = t143 * t146;
	t141 = sin(pkin(6));
	t153 = t145 * t141;
	t152 = t146 * t141;
	t151 = qJD(1) * t145;
	t150 = t142 * qJD(2);
	t144 = cos(pkin(4));
	t149 = t144 * t154;
	t148 = qJD(1) * (-pkin(2) + r_i_i_C(2));
	t147 = t144 * t155 + t152;
	t138 = t147 * qJD(1);
	t136 = -qJD(1) * t149 + t141 * t151;
	t1 = [-(-t149 + t153) * qJD(3) + t146 * t150 + (-t144 * t153 + t154) * t148 + t157 * t138 + (-t158 * t145 * t142 - t146 * pkin(1)) * qJD(1), qJD(1) * t156, -t136, 0; t147 * qJD(3) + t145 * t150 + (t144 * t152 + t155) * t148 + t157 * t136 + (-t145 * pkin(1) + t158 * t156) * qJD(1), t142 * t151, t138, 0; 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:39:02
	% EndTime: 2019-10-09 20:39:02
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (40->24), mult. (126->36), div. (0->0), fcn. (112->6), ass. (0->25)
	t141 = sin(pkin(4));
	t160 = t141 * (pkin(3) + r_i_i_C(1) + qJ(2));
	t159 = -r_i_i_C(2) - qJ(3);
	t140 = sin(pkin(6));
	t144 = sin(qJ(1));
	t158 = t144 * t140;
	t142 = cos(pkin(6));
	t157 = t144 * t142;
	t145 = cos(qJ(1));
	t156 = t145 * t140;
	t155 = t145 * t142;
	t154 = qJD(1) * t144;
	t153 = qJD(1) * t145;
	t152 = t141 * qJD(2);
	t151 = -pkin(2) - r_i_i_C(3) - qJ(4);
	t149 = t140 * t154;
	t148 = t142 * t153;
	t143 = cos(pkin(4));
	t147 = t143 * t157 + t156;
	t146 = t143 * t156 + t157;
	t137 = -t143 * t149 + t148;
	t136 = t147 * qJD(1);
	t135 = t146 * qJD(1);
	t134 = -t143 * t148 + t149;
	t1 = [-t146 * qJD(4) - (-t143 * t155 + t158) * qJD(3) + t145 * t152 + t159 * t136 + t151 * t137 + (-t145 * pkin(1) - t144 * t160) * qJD(1), t141 * t153, -t134, -t135; -(t143 * t158 - t155) * qJD(4) + t147 * qJD(3) + t144 * t152 + t159 * t134 + t151 * t135 + (-t144 * pkin(1) + t145 * t160) * qJD(1), t141 * t154, t136, t137; 0, 0, 0, 0;];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,4);
end