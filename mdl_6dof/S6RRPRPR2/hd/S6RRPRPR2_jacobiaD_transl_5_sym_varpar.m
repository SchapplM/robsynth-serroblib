% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:17
% EndTime: 2019-02-26 21:38:17
% DurationCPUTime: 0.17s
% Computational Cost: add. (270->43), mult. (235->50), div. (0->0), fcn. (161->8), ass. (0->34)
t199 = qJ(2) + pkin(10);
t196 = qJ(4) + t199;
t193 = cos(t196);
t226 = r_i_i_C(3) + qJ(5);
t212 = t226 * t193;
t186 = -sin(qJ(2)) * pkin(2) - pkin(3) * sin(t199);
t176 = t186 * qJD(2);
t192 = sin(t196);
t191 = t192 * qJD(5);
t198 = qJD(2) + qJD(4);
t228 = pkin(4) - r_i_i_C(2);
t233 = (-t228 * t192 + t212) * t198 + (r_i_i_C(1) + pkin(8) + qJ(3) + pkin(7)) * qJD(1) + t176 + t191;
t201 = sin(qJ(1));
t223 = t198 * t201;
t219 = t193 * t223;
t203 = cos(qJ(1));
t221 = qJD(1) * t203;
t232 = t192 * t221 + t219;
t225 = t193 * t198;
t224 = t198 * t192;
t222 = qJD(1) * t201;
t220 = qJD(5) * t193;
t218 = t203 * t225;
t215 = t192 * t222;
t217 = pkin(4) * t215 + r_i_i_C(2) * t218 + t203 * t220;
t213 = t226 * t192;
t210 = -cos(qJ(2)) * pkin(2) - pkin(3) * cos(t199);
t209 = t232 * r_i_i_C(2) + t201 * t220 + t221 * t212;
t208 = -r_i_i_C(2) * t192 - t212;
t207 = -t228 * t224 + t226 * t225 + t191;
t206 = (-pkin(4) * t193 - t213) * t198;
t205 = t210 * qJD(2) + t206;
t204 = qJD(3) + (-t228 * t193 - pkin(1) + t210 - t213) * qJD(1);
t1 = [-t233 * t201 + t204 * t203, t205 * t203 + (-t186 + t208) * t222 + t217, t221, t203 * t206 + t208 * t222 + t217, -t215 + t218, 0; t204 * t201 + t233 * t203 (-pkin(4) * t192 + t186) * t221 + t205 * t201 + t209, t222, -pkin(4) * t219 + (-pkin(4) * t221 - t226 * t223) * t192 + t209, t232, 0; 0, t176 + t207, 0, t207, t224, 0;];
JaD_transl  = t1;
