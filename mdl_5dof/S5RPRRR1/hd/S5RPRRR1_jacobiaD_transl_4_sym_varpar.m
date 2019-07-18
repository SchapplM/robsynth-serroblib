% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RPRRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RPRRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:26:26
% EndTime: 2019-07-18 13:26:26
% DurationCPUTime: 0.10s
% Computational Cost: add. (65->35), mult. (218->71), div. (0->0), fcn. (179->6), ass. (0->27)
t200 = cos(qJ(4));
t202 = cos(qJ(1));
t218 = t200 * t202;
t199 = sin(qJ(1));
t217 = qJD(1) * t199;
t216 = qJD(1) * t202;
t198 = sin(qJ(3));
t215 = qJD(3) * t198;
t201 = cos(qJ(3));
t214 = qJD(3) * t201;
t213 = qJD(3) * t202;
t212 = qJD(4) * t198;
t211 = qJD(4) * t201;
t210 = -qJD(1) + t211;
t209 = qJD(1) * t201 - qJD(4);
t197 = sin(qJ(4));
t208 = r_i_i_C(1) * t200 - r_i_i_C(2) * t197;
t207 = r_i_i_C(1) * t197 + r_i_i_C(2) * t200;
t206 = t210 * t197;
t205 = qJD(3) * t208;
t204 = -r_i_i_C(3) * qJD(3) + t207 * qJD(4);
t203 = t198 * t213 + t209 * t199;
t196 = -t209 * t218 + (t200 * t215 + t206) * t199;
t195 = t210 * t200 * t199 + (-t199 * t215 + t209 * t202) * t197;
t194 = t203 * t200 + t202 * t206;
t193 = t203 * t197 - t210 * t218;
t1 = [-qJ(2) * t217 + t196 * r_i_i_C(1) + t195 * r_i_i_C(2) + t202 * qJD(2) + (-t198 * t216 - t199 * t214) * r_i_i_C(3), t216, (-r_i_i_C(3) * t217 - t202 * t205) * t201 + (t204 * t202 + t208 * t217) * t198, t193 * r_i_i_C(1) + t194 * r_i_i_C(2), 0; qJ(2) * t216 - t194 * r_i_i_C(1) + t193 * r_i_i_C(2) + t199 * qJD(2) + (-t198 * t217 + t201 * t213) * r_i_i_C(3), t217, (r_i_i_C(3) * t216 - t199 * t205) * t201 + (t204 * t199 - t208 * t216) * t198, -t195 * r_i_i_C(1) + t196 * r_i_i_C(2), 0; 0, 0, -t207 * t211 + (r_i_i_C(3) * t201 - t208 * t198) * qJD(3), (t197 * t212 - t200 * t214) * r_i_i_C(2) + (-t197 * t214 - t200 * t212) * r_i_i_C(1), 0;];
JaD_transl  = t1;
