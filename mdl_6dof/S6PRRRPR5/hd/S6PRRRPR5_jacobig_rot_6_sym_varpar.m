% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRRPR5_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_jacobig_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:12:55
% EndTime: 2019-02-26 20:12:55
% DurationCPUTime: 0.06s
% Computational Cost: add. (40->22), mult. (100->47), div. (0->0), fcn. (145->12), ass. (0->30)
t172 = sin(pkin(12));
t174 = sin(pkin(6));
t190 = t172 * t174;
t173 = sin(pkin(7));
t189 = t173 * t174;
t177 = cos(pkin(6));
t188 = t173 * t177;
t175 = cos(pkin(12));
t187 = t175 * t174;
t176 = cos(pkin(7));
t181 = cos(qJ(2));
t186 = t176 * t181;
t179 = sin(qJ(2));
t185 = t177 * t179;
t184 = t177 * t181;
t165 = -t172 * t179 + t175 * t184;
t183 = -t165 * t176 + t173 * t187;
t167 = -t172 * t184 - t175 * t179;
t182 = t167 * t176 + t172 * t189;
t180 = cos(qJ(3));
t178 = sin(qJ(3));
t171 = qJ(4) + pkin(13);
t170 = cos(t171);
t169 = sin(t171);
t168 = -t172 * t185 + t175 * t181;
t166 = t172 * t181 + t175 * t185;
t164 = t177 * t176 - t181 * t189;
t163 = -t167 * t173 + t176 * t190;
t162 = -t165 * t173 - t176 * t187;
t1 = [0, t190, t163, t168 * t178 - t182 * t180, 0 (t168 * t180 + t182 * t178) * t169 - t163 * t170; 0, -t187, t162, t166 * t178 + t183 * t180, 0 (t166 * t180 - t183 * t178) * t169 - t162 * t170; 0, t177, t164, -t180 * t188 + (t178 * t179 - t180 * t186) * t174, 0 (t178 * t188 + (t178 * t186 + t179 * t180) * t174) * t169 - t164 * t170;];
Jg_rot  = t1;
