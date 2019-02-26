% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR1_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobigD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:37
% EndTime: 2019-02-26 20:18:37
% DurationCPUTime: 0.04s
% Computational Cost: add. (46->12), mult. (72->30), div. (0->0), fcn. (74->8), ass. (0->21)
t178 = qJ(3) + qJ(4) + qJ(5);
t175 = sin(t178);
t180 = sin(pkin(6));
t191 = t180 * t175;
t182 = cos(pkin(6));
t183 = sin(qJ(2));
t190 = t182 * t183;
t184 = cos(qJ(2));
t189 = t182 * t184;
t188 = qJD(2) * t175;
t187 = qJD(2) * t180;
t179 = sin(pkin(12));
t181 = cos(pkin(12));
t186 = t179 * t184 + t181 * t190;
t185 = -t179 * t190 + t181 * t184;
t177 = qJD(3) + qJD(4) + qJD(5);
t176 = cos(t178);
t174 = t183 * t187;
t173 = t185 * qJD(2);
t172 = t186 * qJD(2);
t1 = [0, 0, t173, t173, t173 (t185 * t176 + t179 * t191) * t177 + (-t179 * t189 - t181 * t183) * t188; 0, 0, t172, t172, t172 (t186 * t176 - t181 * t191) * t177 + (-t179 * t183 + t181 * t189) * t188; 0, 0, t174, t174, t174, t180 * t183 * t177 * t176 + (t177 * t182 + t184 * t187) * t175;];
JgD_rot  = t1;
