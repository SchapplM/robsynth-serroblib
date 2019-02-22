% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:25
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PPRRRP2_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_jacobiR_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:25:25
% EndTime: 2019-02-22 09:25:25
% DurationCPUTime: 0.19s
% Computational Cost: add. (175->43), mult. (516->91), div. (0->0), fcn. (712->14), ass. (0->46)
t206 = cos(qJ(3));
t180 = sin(pkin(11));
t186 = cos(pkin(6));
t205 = t180 * t186;
t181 = sin(pkin(7));
t204 = t181 * t186;
t182 = sin(pkin(6));
t203 = t182 * t181;
t185 = cos(pkin(7));
t202 = t182 * t185;
t183 = cos(pkin(12));
t201 = t183 * t185;
t184 = cos(pkin(11));
t200 = t184 * t186;
t187 = sin(qJ(5));
t191 = cos(qJ(4));
t199 = t187 * t191;
t190 = cos(qJ(5));
t198 = t190 * t191;
t197 = t182 * t206;
t196 = t181 * t197;
t179 = sin(pkin(12));
t195 = -t180 * t179 + t183 * t200;
t194 = t184 * t179 + t183 * t205;
t193 = t195 * t185;
t192 = t194 * t185;
t189 = sin(qJ(3));
t188 = sin(qJ(4));
t175 = -t179 * t205 + t184 * t183;
t174 = t179 * t200 + t180 * t183;
t173 = -t183 * t203 + t186 * t185;
t170 = t180 * t202 + t194 * t181;
t169 = -t195 * t181 - t184 * t202;
t168 = t189 * t204 + (t206 * t179 + t189 * t201) * t182;
t167 = t182 * t179 * t189 - t197 * t201 - t206 * t204;
t166 = t168 * t191 + t173 * t188;
t165 = -t168 * t188 + t173 * t191;
t164 = t175 * t206 + (t180 * t203 - t192) * t189;
t163 = t175 * t189 - t180 * t196 + t206 * t192;
t162 = t174 * t206 + (-t184 * t203 + t193) * t189;
t161 = t174 * t189 + t184 * t196 - t206 * t193;
t160 = t164 * t191 + t170 * t188;
t159 = -t164 * t188 + t170 * t191;
t158 = t162 * t191 + t169 * t188;
t157 = -t162 * t188 + t169 * t191;
t1 = [0, 0, -t163 * t198 + t164 * t187, t159 * t190, -t160 * t187 + t163 * t190, 0; 0, 0, -t161 * t198 + t162 * t187, t157 * t190, -t158 * t187 + t161 * t190, 0; 0, 0, -t167 * t198 + t168 * t187, t165 * t190, -t166 * t187 + t167 * t190, 0; 0, 0, t163 * t199 + t164 * t190, -t159 * t187, -t160 * t190 - t163 * t187, 0; 0, 0, t161 * t199 + t162 * t190, -t157 * t187, -t158 * t190 - t161 * t187, 0; 0, 0, t167 * t199 + t168 * t190, -t165 * t187, -t166 * t190 - t167 * t187, 0; 0, 0, -t163 * t188, t160, 0, 0; 0, 0, -t161 * t188, t158, 0, 0; 0, 0, -t167 * t188, t166, 0, 0;];
JR_rot  = t1;
