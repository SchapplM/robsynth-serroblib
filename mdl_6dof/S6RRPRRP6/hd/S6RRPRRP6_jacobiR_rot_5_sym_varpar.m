% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:34
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRRP6_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:34:01
% EndTime: 2019-02-22 11:34:01
% DurationCPUTime: 0.12s
% Computational Cost: add. (154->31), mult. (425->65), div. (0->0), fcn. (606->12), ass. (0->37)
t173 = sin(qJ(4));
t177 = cos(qJ(4));
t171 = cos(pkin(6));
t168 = sin(pkin(11));
t170 = cos(pkin(11));
t174 = sin(qJ(2));
t178 = cos(qJ(2));
t181 = t178 * t168 + t174 * t170;
t162 = t181 * t171;
t163 = t174 * t168 - t178 * t170;
t175 = sin(qJ(1));
t179 = cos(qJ(1));
t183 = t179 * t162 - t175 * t163;
t169 = sin(pkin(6));
t188 = t169 * t179;
t148 = t173 * t188 - t177 * t183;
t180 = t163 * t171;
t152 = -t175 * t181 - t179 * t180;
t172 = sin(qJ(5));
t176 = cos(qJ(5));
t193 = t148 * t172 - t152 * t176;
t192 = t148 * t176 + t152 * t172;
t189 = t169 * t175;
t187 = t172 * t177;
t185 = t176 * t177;
t182 = -t175 * t162 - t179 * t163;
t146 = -t173 * t183 - t177 * t188;
t161 = t181 * t169;
t160 = t163 * t169;
t158 = t161 * t177 + t171 * t173;
t157 = -t161 * t173 + t171 * t177;
t155 = t175 * t180 - t179 * t181;
t150 = t173 * t189 + t177 * t182;
t149 = t173 * t182 - t177 * t189;
t145 = t150 * t176 - t155 * t172;
t144 = -t150 * t172 - t155 * t176;
t1 = [t192, t155 * t185 + t172 * t182, 0, -t149 * t176, t144, 0; t145, t152 * t185 + t172 * t183, 0, t146 * t176, t193, 0; 0, -t160 * t185 + t161 * t172, 0, t157 * t176, -t158 * t172 + t160 * t176, 0; -t193, -t155 * t187 + t176 * t182, 0, t149 * t172, -t145, 0; t144, -t152 * t187 + t176 * t183, 0, -t146 * t172, t192, 0; 0, t160 * t187 + t161 * t176, 0, -t157 * t172, -t158 * t176 - t160 * t172, 0; t146, t155 * t173, 0, t150, 0, 0; t149, t152 * t173, 0, -t148, 0, 0; 0, -t160 * t173, 0, t158, 0, 0;];
JR_rot  = t1;
