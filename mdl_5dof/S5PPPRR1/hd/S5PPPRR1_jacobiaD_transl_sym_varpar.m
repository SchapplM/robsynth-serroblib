% Zeitableitung der analytischen Jacobi-Matrix (Translatorisch) für beliebiges Segment von
% S5PPPRR1
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% JaD_transl [3x5]
%   Translatorischer Teil der analytischen Jacobi-Matrix (Zeitableitung)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:17
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5PPPRR1_jacobiaD_transl_sym_varpar(qJ, qJD, link_index, r_i_i_C, ...
  pkin)


%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_jacobiaD_transl_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_jacobiaD_transl_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5PPPRR1_jacobiaD_transl_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PPPRR1_jacobiaD_transl_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_jacobiaD_transl_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_transl_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_transl_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_transl_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_transl_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_transl = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_transl_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:16
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (17->8), mult. (30->20), div. (0->0), fcn. (24->6), ass. (0->9)
	t61 = sin(pkin(7));
	t62 = cos(pkin(8));
	t65 = t61 * t62;
	t63 = cos(pkin(7));
	t64 = t62 * t63;
	t59 = pkin(9) + qJ(4);
	t58 = cos(t59);
	t57 = sin(t59);
	t1 = [0, 0, 0, ((-t57 * t61 - t58 * t64) * r_i_i_C(1) + (t57 * t64 - t58 * t61) * r_i_i_C(2)) * qJD(4), 0; 0, 0, 0, ((t57 * t63 - t58 * t65) * r_i_i_C(1) + (t57 * t65 + t58 * t63) * r_i_i_C(2)) * qJD(4), 0; 0, 0, 0, (-r_i_i_C(1) * t58 + r_i_i_C(2) * t57) * sin(pkin(8)) * qJD(4), 0;];
	JaD_transl = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_transl_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:17:16
	% EndTime: 2019-10-24 10:17:17
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (110->27), mult. (191->56), div. (0->0), fcn. (172->8), ass. (0->27)
	t186 = pkin(9) + qJ(4);
	t185 = cos(t186);
	t188 = sin(pkin(7));
	t189 = cos(pkin(8));
	t184 = sin(t186);
	t190 = cos(pkin(7));
	t201 = t190 * t184;
	t207 = -t188 * t185 + t189 * t201;
	t206 = -pkin(6) - r_i_i_C(3);
	t187 = sin(pkin(8));
	t191 = sin(qJ(5));
	t205 = t187 * t191;
	t192 = cos(qJ(5));
	t204 = t187 * t192;
	t202 = t188 * t189;
	t200 = t190 * t185;
	t198 = t191 * r_i_i_C(1) + t192 * r_i_i_C(2);
	t197 = t192 * r_i_i_C(1) - t191 * r_i_i_C(2) + pkin(4);
	t196 = t198 * t184;
	t182 = t188 * t184 + t189 * t200;
	t180 = t185 * t202 - t201;
	t195 = t184 * t202 + t200;
	t194 = qJD(5) * t198;
	t193 = qJD(4) * t197;
	t177 = t207 * qJD(4);
	t175 = t195 * qJD(4);
	t1 = [0, 0, 0, t206 * t177 - t182 * t193 + t207 * t194, t198 * t177 + ((-t182 * t192 - t190 * t205) * r_i_i_C(1) + (t182 * t191 - t190 * t204) * r_i_i_C(2)) * qJD(5); 0, 0, 0, t206 * t175 - t180 * t193 + t195 * t194, t198 * t175 + ((-t180 * t192 - t188 * t205) * r_i_i_C(1) + (t180 * t191 - t188 * t204) * r_i_i_C(2)) * qJD(5); 0, 0, 0, (qJD(5) * t196 + (t206 * t184 - t197 * t185) * qJD(4)) * t187, t187 * qJD(4) * t196 + ((-t185 * t204 + t189 * t191) * r_i_i_C(1) + (t185 * t205 + t189 * t192) * r_i_i_C(2)) * qJD(5);];
	JaD_transl = t1;
else
	JaD_transl=NaN(3,5);
end