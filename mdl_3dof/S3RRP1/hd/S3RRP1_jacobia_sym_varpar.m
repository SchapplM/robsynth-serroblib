% Analytische Jacobi-Matrix für beliebiges Segment von
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S3RRP1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% 
% Output:
% Ja [6x3]
%   Analytische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:12
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja = S3RRP1_jacobia_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),uint8(0),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_jacobia_sym_varpar: qJ has to be [3x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S3RRP1_jacobia_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S3RRP1_jacobia_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_jacobia_sym_varpar: pkin has to be [4x1] (double)');

% Function calls
Ja_transl = S3RRP1_jacobia_transl_sym_varpar(qJ, link_index, r_i_i_C, ...
  pkin);
Ja_rot = S3RRP1_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin);

Ja = [Ja_transl; Ja_rot];
end