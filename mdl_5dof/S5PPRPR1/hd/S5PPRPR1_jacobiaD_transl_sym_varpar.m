% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPRPR1
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% (Ist für translatorischen Teil egal, kennzeichnet nur den Rechenweg der Herleitung)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt (0=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:18
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPRPR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPRPR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPRPR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->3), mult. (16->8), div. (0->0), fcn. (10->4), ass. (0->5)
	t14 = pkin(8) + qJ(3);
	t12 = sin(t14);
	t13 = cos(t14);
	t17 = qJD(3) * (-r_i_i_C(1) * t13 + r_i_i_C(2) * t12);
	t1 = [0, 0, cos(pkin(7)) * t17, 0, 0; 0, 0, sin(pkin(7)) * t17, 0, 0; 0, 0, (-r_i_i_C(1) * t12 - r_i_i_C(2) * t13) * qJD(3), 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (36->8), mult. (56->16), div. (0->0), fcn. (41->6), ass. (0->10)
	t128 = r_i_i_C(3) + qJ(4);
	t120 = pkin(8) + qJ(3);
	t119 = cos(t120);
	t127 = qJD(3) * t119;
	t126 = -r_i_i_C(1) * cos(pkin(9)) + r_i_i_C(2) * sin(pkin(9)) - pkin(3);
	t118 = sin(t120);
	t125 = qJD(4) * t119 + (-t128 * t118 + t126 * t119) * qJD(3);
	t124 = cos(pkin(7));
	t122 = sin(pkin(7));
	t1 = [0, 0, t125 * t124, t124 * t127, 0; 0, 0, t125 * t122, t122 * t127, 0; 0, 0, t118 * qJD(4) + (t126 * t118 + t128 * t119) * qJD(3), qJD(3) * t118, 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:18:05
	% EndTime: 2019-10-24 10:18:05
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (107->24), mult. (129->49), div. (0->0), fcn. (100->7), ass. (0->19)
	t171 = r_i_i_C(3) + pkin(6) + qJ(4);
	t158 = pkin(8) + qJ(3);
	t156 = cos(t158);
	t159 = sin(pkin(7));
	t170 = t156 * t159;
	t160 = cos(pkin(7));
	t169 = t156 * t160;
	t168 = qJD(3) * t156;
	t154 = sin(t158);
	t167 = qJD(5) * t154;
	t157 = pkin(9) + qJ(5);
	t153 = sin(t157);
	t155 = cos(t157);
	t166 = r_i_i_C(1) * t153 + r_i_i_C(2) * t155;
	t165 = -r_i_i_C(1) * t155 + r_i_i_C(2) * t153 - cos(pkin(9)) * pkin(4) - pkin(3);
	t164 = t166 * t154;
	t163 = qJD(3) * t164;
	t162 = qJD(4) * t156 + qJD(5) * t164 + (-t154 * t171 + t165 * t156) * qJD(3);
	t1 = [0, 0, t162 * t160, t160 * t168, t160 * t163 + ((-t153 * t159 - t155 * t169) * r_i_i_C(1) + (t153 * t169 - t155 * t159) * r_i_i_C(2)) * qJD(5); 0, 0, t162 * t159, t159 * t168, t159 * t163 + ((t153 * t160 - t155 * t170) * r_i_i_C(1) + (t153 * t170 + t155 * t160) * r_i_i_C(2)) * qJD(5); 0, 0, t154 * qJD(4) - t166 * t156 * qJD(5) + (t165 * t154 + t156 * t171) * qJD(3), qJD(3) * t154, (t153 * t167 - t155 * t168) * r_i_i_C(2) + (-t153 * t168 - t155 * t167) * r_i_i_C(1);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end