% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 11:26
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRPRPR7_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR7_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR7_jacobiR_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 11:26:19
% EndTime: 2019-02-22 11:26:19
% DurationCPUTime: 0.04s
% Computational Cost: add. (50->13), mult. (56->8), div. (0->0), fcn. (90->6), ass. (0->14)
t43 = qJ(4) + pkin(10);
t41 = sin(t43);
t42 = cos(t43);
t44 = sin(qJ(2));
t46 = cos(qJ(2));
t49 = t46 * t41 - t44 * t42;
t48 = t44 * t41 + t46 * t42;
t47 = cos(qJ(1));
t45 = sin(qJ(1));
t37 = t48 * t47;
t36 = t49 * t47;
t35 = t48 * t45;
t34 = t49 * t45;
t1 = [-t35, t36, 0, -t36, 0, 0; t37, t34, 0, -t34, 0, 0; 0, t48, 0, -t48, 0, 0; t34, t37, 0, -t37, 0, 0; -t36, t35, 0, -t35, 0, 0; 0, -t49, 0, t49, 0, 0; -t47, 0, 0, 0, 0, 0; -t45, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
