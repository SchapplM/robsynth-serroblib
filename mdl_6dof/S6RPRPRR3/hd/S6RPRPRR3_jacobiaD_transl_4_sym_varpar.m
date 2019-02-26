% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:50:20
% EndTime: 2019-02-26 20:50:20
% DurationCPUTime: 0.13s
% Computational Cost: add. (113->22), mult. (182->34), div. (0->0), fcn. (135->8), ass. (0->19)
t167 = sin(qJ(3));
t168 = cos(qJ(3));
t165 = sin(pkin(11));
t166 = cos(pkin(11));
t176 = r_i_i_C(1) * t166 - r_i_i_C(2) * t165 + pkin(3);
t180 = r_i_i_C(3) + qJ(4);
t173 = t176 * t167 - t180 * t168;
t182 = qJD(1) * t173;
t181 = t173 * qJD(3) - t167 * qJD(4);
t179 = qJD(1) * t167;
t178 = qJD(3) * t168;
t175 = t165 * r_i_i_C(1) + t166 * r_i_i_C(2) + pkin(7);
t174 = -t180 * t167 - t176 * t168;
t171 = -pkin(2) + t174;
t170 = t174 * qJD(3) + qJD(4) * t168;
t164 = qJ(1) + pkin(10);
t163 = cos(t164);
t162 = sin(t164);
t1 = [t181 * t162 + (-cos(qJ(1)) * pkin(1) - t175 * t162 + t171 * t163) * qJD(1), 0, t162 * t182 + t170 * t163, -t162 * t179 + t163 * t178, 0, 0; -t181 * t163 + (-sin(qJ(1)) * pkin(1) + t175 * t163 + t171 * t162) * qJD(1), 0, t170 * t162 - t163 * t182, t162 * t178 + t163 * t179, 0, 0; 0, 0, -t181, qJD(3) * t167, 0, 0;];
JaD_transl  = t1;
