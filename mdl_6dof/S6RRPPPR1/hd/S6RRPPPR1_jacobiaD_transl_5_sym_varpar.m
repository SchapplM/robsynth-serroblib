% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:21:59
% EndTime: 2019-02-26 21:21:59
% DurationCPUTime: 0.16s
% Computational Cost: add. (184->41), mult. (312->62), div. (0->0), fcn. (248->8), ass. (0->32)
t189 = qJ(2) + pkin(9);
t187 = sin(t189);
t188 = cos(t189);
t190 = sin(pkin(10));
t208 = qJD(5) * t190;
t202 = t187 * qJD(4) + t188 * t208;
t218 = r_i_i_C(2) + qJ(4);
t204 = t218 * t188 - sin(qJ(2)) * pkin(2);
t227 = (-pkin(3) * t187 + t204) * qJD(2) - qJD(1) * (-qJ(3) - pkin(7)) + t202;
t191 = cos(pkin(10));
t217 = r_i_i_C(3) + qJ(5);
t222 = pkin(4) + r_i_i_C(1);
t223 = t217 * t190 + t222 * t191 + pkin(3);
t225 = t223 * t187 - t204;
t224 = -t218 * t187 - cos(qJ(2)) * pkin(2);
t194 = sin(qJ(1));
t216 = t190 * t194;
t196 = cos(qJ(1));
t215 = t190 * t196;
t214 = t194 * t191;
t213 = t196 * t191;
t212 = qJD(1) * t194;
t211 = qJD(1) * t196;
t210 = qJD(2) * t194;
t209 = qJD(2) * t196;
t207 = t187 * t210;
t206 = t187 * t209;
t199 = -t191 * qJD(5) + qJD(3) + (-pkin(3) * t188 - pkin(1) + t224) * qJD(1);
t197 = -t187 * t208 + qJD(4) * t188 + (-t188 * t223 + t224) * qJD(2);
t184 = -t190 * t207 + (t188 * t215 - t214) * qJD(1);
t182 = t190 * t206 + (t188 * t216 + t213) * qJD(1);
t1 = [t222 * (t191 * t207 + (-t188 * t213 - t216) * qJD(1)) - t217 * t184 + t199 * t196 - t227 * t194, t197 * t196 + t225 * t212, t211, -t187 * t212 + t188 * t209, -t182, 0; t222 * (-t191 * t206 + (-t188 * t214 + t215) * qJD(1)) - t217 * t182 + t199 * t194 + t227 * t196, t197 * t194 - t211 * t225, t212, t187 * t211 + t188 * t210, t184, 0; 0, -qJD(2) * t225 + t202, 0, qJD(2) * t187, qJD(2) * t188 * t190, 0;];
JaD_transl  = t1;
