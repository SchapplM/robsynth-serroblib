% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR4_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR4_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:05:53
% EndTime: 2019-02-26 20:05:53
% DurationCPUTime: 0.17s
% Computational Cost: add. (124->35), mult. (404->66), div. (0->0), fcn. (384->8), ass. (0->32)
t231 = sin(pkin(11));
t233 = cos(pkin(11));
t236 = sin(qJ(2));
t234 = cos(pkin(6));
t238 = cos(qJ(2));
t248 = t234 * t238;
t259 = -t231 * t236 + t233 * t248;
t249 = t234 * t236;
t226 = t231 * t238 + t233 * t249;
t237 = cos(qJ(3));
t232 = sin(pkin(6));
t235 = sin(qJ(3));
t251 = t232 * t235;
t258 = -t226 * t237 + t233 * t251;
t254 = r_i_i_C(3) + qJ(4);
t256 = pkin(3) + r_i_i_C(1);
t257 = t254 * t235 + t256 * t237 + pkin(2);
t255 = pkin(8) + r_i_i_C(2);
t250 = t232 * t237;
t245 = qJD(2) * t232 * t238;
t242 = t231 * t249 - t233 * t238;
t244 = t231 * t251 - t237 * t242;
t243 = t231 * t248 + t233 * t236;
t241 = t234 * t235 + t236 * t250;
t240 = qJD(2) * t257;
t239 = qJD(4) * t235 + (-t256 * t235 + t254 * t237) * qJD(3);
t223 = t243 * qJD(2);
t221 = t259 * qJD(2);
t219 = t241 * qJD(3) + t235 * t245;
t217 = t244 * qJD(3) - t223 * t235;
t215 = -t258 * qJD(3) + t221 * t235;
t1 = [0, -t255 * t223 - t239 * t243 + t242 * t240, t244 * qJD(4) + t254 * (-t223 * t237 + (t231 * t250 + t235 * t242) * qJD(3)) - t256 * t217, t217, 0, 0; 0, t255 * t221 - t226 * t240 + t239 * t259, -t258 * qJD(4) + t254 * (t221 * t237 + (-t226 * t235 - t233 * t250) * qJD(3)) - t256 * t215, t215, 0, 0; 0 (t239 * t238 + (-t257 * t236 + t255 * t238) * qJD(2)) * t232, t241 * qJD(4) + t254 * (t237 * t245 + (t234 * t237 - t236 * t251) * qJD(3)) - t256 * t219, t219, 0, 0;];
JaD_transl  = t1;
