% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:34
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRPRRP2_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:34:04
% EndTime: 2019-02-22 09:34:04
% DurationCPUTime: 0.14s
% Computational Cost: add. (111->25), mult. (319->65), div. (0->0), fcn. (454->12), ass. (0->33)
t169 = sin(pkin(11));
t172 = cos(pkin(11));
t177 = sin(qJ(2));
t180 = cos(qJ(2));
t165 = t177 * t169 - t180 * t172;
t171 = sin(pkin(6));
t176 = sin(qJ(4));
t190 = t171 * t176;
t179 = cos(qJ(4));
t189 = t171 * t179;
t175 = sin(qJ(5));
t188 = t175 * t179;
t178 = cos(qJ(5));
t186 = t178 * t179;
t174 = cos(pkin(6));
t182 = t180 * t169 + t177 * t172;
t164 = t182 * t174;
t170 = sin(pkin(10));
t173 = cos(pkin(10));
t184 = t173 * t164 - t170 * t165;
t183 = -t170 * t164 - t173 * t165;
t181 = t165 * t174;
t163 = t182 * t171;
t162 = t165 * t171;
t160 = t163 * t179 + t174 * t176;
t159 = -t163 * t176 + t174 * t179;
t157 = t170 * t181 - t173 * t182;
t154 = -t170 * t182 - t173 * t181;
t152 = t170 * t190 + t179 * t183;
t151 = t170 * t189 - t176 * t183;
t150 = -t173 * t190 + t179 * t184;
t149 = -t173 * t189 - t176 * t184;
t1 = [0, t157 * t186 + t175 * t183, 0, t151 * t178, -t152 * t175 - t157 * t178, 0; 0, t154 * t186 + t175 * t184, 0, t149 * t178, -t150 * t175 - t154 * t178, 0; 0, -t162 * t186 + t163 * t175, 0, t159 * t178, -t160 * t175 + t162 * t178, 0; 0, t157 * t176, 0, t152, 0, 0; 0, t154 * t176, 0, t150, 0, 0; 0, -t162 * t176, 0, t160, 0, 0; 0, t157 * t188 - t178 * t183, 0, t151 * t175, t152 * t178 - t157 * t175, 0; 0, t154 * t188 - t178 * t184, 0, t149 * t175, t150 * t178 - t154 * t175, 0; 0, -t162 * t188 - t163 * t178, 0, t159 * t175, t160 * t178 + t162 * t175, 0;];
JR_rot  = t1;
