% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:22:36
% EndTime: 2019-02-26 21:22:37
% DurationCPUTime: 0.14s
% Computational Cost: add. (160->27), mult. (252->40), div. (0->0), fcn. (189->8), ass. (0->20)
t174 = sin(pkin(10));
t175 = cos(pkin(10));
t173 = qJ(2) + pkin(9);
t171 = sin(t173);
t172 = cos(t173);
t188 = r_i_i_C(1) * t174 + r_i_i_C(2) * t175 + qJ(4);
t190 = pkin(3) + r_i_i_C(3) + qJ(5);
t184 = t190 * t171 - t188 * t172 + sin(qJ(2)) * pkin(2);
t181 = -t184 * qJD(2) + t171 * qJD(4) + t172 * qJD(5);
t199 = t181 + (t175 * r_i_i_C(1) - t174 * r_i_i_C(2) + pkin(4) + pkin(7) + qJ(3)) * qJD(1);
t198 = -t188 * t171 - t190 * t172 - cos(qJ(2)) * pkin(2);
t178 = sin(qJ(1));
t194 = qJD(1) * t178;
t180 = cos(qJ(1));
t193 = qJD(1) * t180;
t192 = qJD(2) * t178;
t191 = qJD(2) * t180;
t183 = qJD(3) + (-pkin(1) + t198) * qJD(1);
t182 = t198 * qJD(2) + qJD(4) * t172 - qJD(5) * t171;
t1 = [-t199 * t178 + t183 * t180, t182 * t180 + t184 * t194, t193, -t171 * t194 + t172 * t191, -t171 * t191 - t172 * t194, 0; t183 * t178 + t199 * t180, t182 * t178 - t184 * t193, t194, t171 * t193 + t172 * t192, -t171 * t192 + t172 * t193, 0; 0, t181, 0, qJD(2) * t171, qJD(2) * t172, 0;];
JaD_transl  = t1;
