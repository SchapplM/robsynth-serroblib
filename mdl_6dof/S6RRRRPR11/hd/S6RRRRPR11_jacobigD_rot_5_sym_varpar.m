% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRPR11_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:22
% EndTime: 2019-02-26 22:36:22
% DurationCPUTime: 0.05s
% Computational Cost: add. (22->16), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->21)
t159 = sin(pkin(6));
t164 = cos(qJ(3));
t178 = t159 * t164;
t166 = cos(qJ(1));
t177 = t159 * t166;
t162 = sin(qJ(2));
t163 = sin(qJ(1));
t176 = t162 * t163;
t175 = t162 * t166;
t165 = cos(qJ(2));
t174 = t163 * t165;
t173 = t166 * t165;
t172 = qJD(1) * t159;
t161 = sin(qJ(3));
t171 = qJD(2) * t161;
t160 = cos(pkin(6));
t170 = t160 * t173 - t176;
t169 = t160 * t174 + t175;
t168 = t160 * t175 + t174;
t167 = -t160 * t176 + t173;
t1 = [0, t166 * t172, t170 * qJD(1) + t167 * qJD(2) (t163 * t159 * t161 + t167 * t164) * qJD(3) - t169 * t171 + (-t168 * t161 - t164 * t177) * qJD(1), 0, 0; 0, t163 * t172, t169 * qJD(1) + t168 * qJD(2) (-t161 * t177 + t168 * t164) * qJD(3) + t170 * t171 + (t167 * t161 - t163 * t178) * qJD(1), 0, 0; 0, 0, t159 * qJD(2) * t162, t159 * t165 * t171 + (t160 * t161 + t162 * t178) * qJD(3), 0, 0;];
JgD_rot  = t1;
