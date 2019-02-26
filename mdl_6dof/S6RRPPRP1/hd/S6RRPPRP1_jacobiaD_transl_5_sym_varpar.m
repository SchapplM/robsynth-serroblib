% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRP1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:46
% EndTime: 2019-02-26 21:24:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (264->46), mult. (329->72), div. (0->0), fcn. (261->10), ass. (0->37)
t215 = cos(pkin(10)) * pkin(4) + pkin(3);
t222 = qJ(2) + pkin(9);
t218 = sin(t222);
t220 = cos(t222);
t253 = r_i_i_C(3) + pkin(8) + qJ(4);
t235 = t253 * t220 - sin(qJ(2)) * pkin(2);
t244 = t218 * qJD(4);
t262 = (-t215 * t218 + t235) * qJD(2) + (pkin(4) * sin(pkin(10)) + qJ(3) + pkin(7)) * qJD(1) + t244;
t221 = pkin(10) + qJ(5);
t217 = sin(t221);
t219 = cos(t221);
t236 = r_i_i_C(1) * t219 - r_i_i_C(2) * t217 + t215;
t231 = -t236 * t218 + t235;
t229 = cos(qJ(1));
t245 = qJD(5) * t220;
t240 = -qJD(1) + t245;
t260 = t229 * t240;
t258 = -t253 * t218 - cos(qJ(2)) * pkin(2);
t239 = qJD(1) * t220 - qJD(5);
t227 = sin(qJ(1));
t248 = qJD(2) * t227;
t257 = -t218 * t248 + t239 * t229;
t251 = qJD(1) * t227;
t250 = qJD(1) * t229;
t249 = qJD(2) * t220;
t247 = qJD(2) * t229;
t246 = qJD(5) * t218;
t238 = r_i_i_C(1) * t217 + r_i_i_C(2) * t219;
t237 = t240 * t227;
t233 = t218 * t247 + t239 * t227;
t232 = qJD(3) + (-t215 * t220 - pkin(1) + t258) * qJD(1);
t230 = qJD(4) * t220 + t238 * t246 + (-t236 * t220 + t258) * qJD(2);
t214 = t217 * t237 - t257 * t219;
t213 = t257 * t217 + t219 * t237;
t212 = t217 * t260 + t233 * t219;
t211 = t233 * t217 - t219 * t260;
t1 = [t214 * r_i_i_C(1) + t213 * r_i_i_C(2) - t262 * t227 + t232 * t229, t230 * t229 - t231 * t251, t250, -t218 * t251 + t220 * t247, t211 * r_i_i_C(1) + t212 * r_i_i_C(2), 0; -t212 * r_i_i_C(1) + t211 * r_i_i_C(2) + t232 * t227 + t262 * t229, t230 * t227 + t231 * t250, t251, t218 * t250 + t220 * t248, -t213 * r_i_i_C(1) + t214 * r_i_i_C(2), 0; 0, t231 * qJD(2) - t238 * t245 + t244, 0, qJD(2) * t218 (t217 * t246 - t219 * t249) * r_i_i_C(2) + (-t217 * t249 - t219 * t246) * r_i_i_C(1), 0;];
JaD_transl  = t1;
